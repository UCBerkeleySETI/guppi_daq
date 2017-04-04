/* guppi_pktsock_thread_codd.c
 *
 * Routine to read packets from network and put them
 * into shared memory blocks.
 *
 * This one is specific to coherent dedispersion mode because
 * it has the capability to overlap the output data blocks
 * and also performs a transpose on the input.
 */

#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <sched.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <endian.h>

#include "fitshead.h"
#include "guppi_params.h"
#include "psrfits.h"
#include "guppi_error.h"
#include "guppi_status.h"
#include "guppi_databuf.h"
#include "guppi_udp.h"
#include "guppi_time.h"

#include "hashpipe.h"

#define STATUS_KEY "NETSTAT"  /* Define before guppi_threads.h */
#include "guppi_threads.h"

#define PKTSOCK_BYTES_PER_FRAME (16384)
#define PKTSOCK_FRAMES_PER_BLOCK (8)
#define PKTSOCK_NBLOCKS (800)
#define PKTSOCK_NFRAMES (PKTSOCK_FRAMES_PER_BLOCK * PKTSOCK_NBLOCKS)

/* It's easier to just make these global ... */
static unsigned long long npacket_total=0, ndropped_total=0, nbogus_total=0;

/* Structs/functions to more easily deal with multiple 
 * active blocks being filled
 */
struct datablock_stats {
    struct guppi_databuf *db;      // Pointer to overall shared mem databuf
    int block_idx;                 // Block index number in databuf
    unsigned long long packet_idx; // Index of first packet number in block
    size_t packet_data_size;       // Data size of each packet
    int packets_per_block;         // Total number of packets to go in the block
    int overlap_packets;           // Overlap between blocks in packets
    int npacket;                   // Number of packets filled so far
    int ndropped;                  // Number of dropped packets so far
    unsigned long long last_pkt;   // Last packet seq number written to block
};

// Defined in guppi_net_thread_codd.c

/* get the thread specific pid */
pid_t gettid();

/* Reset all counters */
void reset_stats(struct datablock_stats *d);

/* Reset block params */
void reset_block(struct datablock_stats *d);

/* Initialize block struct */
void init_block(struct datablock_stats *d, struct guppi_databuf *db, 
        size_t packet_data_size, int packets_per_block, int overlap_packets);

/* Update block header info, set filled status */
void finalize_block(struct datablock_stats *d);

/* Push all blocks down a level, losing the first one */
void block_stack_push(struct datablock_stats *d, int nblock);

/* Go to next block in set */
void increment_block(struct datablock_stats *d, 
        unsigned long long next_seq_num);

/* Check whether a certain seq num belongs in the data block */
int block_packet_check(struct datablock_stats *d, 
        unsigned long long seq_num);

/* Return packet index from a pktsock frame that is assumed to contain a UDP
 * packet.
 */
unsigned long long guppi_pktsock_seq_num(const unsigned char *p_frame) {
    // XXX Temp for new baseband mode, blank out top 8 bits which 
    // contain channel info.
    unsigned long long tmp = be64toh(*(uint64_t *)PKT_UDP_DATA(p_frame));
    tmp &= 0x00FFFFFFFFFFFFFF;
    return tmp ;
}

/* Write a search mode (filterbank) style packet (from pktsock) into the
 * datablock.  Also zeroes out any dropped packets.
 */
void write_search_packet_to_block_from_pktsock_frame(
        struct datablock_stats *d, unsigned char *p_frame) {
    const unsigned long long seq_num = guppi_pktsock_seq_num(p_frame);
    int next_pos = seq_num - d->packet_idx;
    int cur_pos=0;
    if (d->last_pkt > d->packet_idx) cur_pos = d->last_pkt - d->packet_idx + 1;
    char *dataptr = guppi_databuf_data(d->db, d->block_idx) 
        + cur_pos*d->packet_data_size;
    for (; cur_pos<next_pos; cur_pos++) {
        memset(dataptr, 0, d->packet_data_size);
        dataptr += d->packet_data_size;
        d->npacket++;
        d->ndropped++;
    }
    guppi_udp_packet_data_copy_from_payload(dataptr,
        (char *)PKT_UDP_DATA(p_frame), (size_t)PKT_UDP_SIZE(p_frame));
    d->last_pkt = seq_num;
    //d->packet_idx++; // XXX I think this is wrong..
    d->npacket++;
}

/* Write a baseband mode packet into the block.  Includes a 
 * corner-turn (aka transpose) of dimension nchan.
 */
void write_baseband_packet_to_block_from_pktsock_frame(
        struct datablock_stats *d, unsigned char *p_frame, int nchan) {

    const unsigned long long seq_num = guppi_pktsock_seq_num(p_frame);
    int block_pkt_idx = seq_num - d->packet_idx;
    guppi_udp_packet_data_copy_transpose_from_payload(
            guppi_databuf_data(d->db, d->block_idx),
            nchan, block_pkt_idx, d->packets_per_block,
            (char *)PKT_UDP_DATA(p_frame),
            (size_t)PKT_UDP_SIZE(p_frame));

    /* Consider any skipped packets to have been dropped,
     * update counters.
     */
    if (d->last_pkt < d->packet_idx) d->last_pkt = d->packet_idx;

    if (seq_num == d->last_pkt) {
        d->npacket++;
    } else {
        d->npacket += seq_num - d->last_pkt;
        d->ndropped += seq_num - d->last_pkt - 1;
    }

    d->last_pkt = seq_num;
}

/* This thread is passed a single arg, pointer
 * to the guppi_pktsock_params struct.  This thread should 
 * be cancelled and restarted if any hardware params
 * change, as this potentially affects packet size, etc.
 */
void *guppi_pktsock_thread_codd(void *_args) {

    guppi_warn("guppi_pktsock_thread_codd", "Thread starting!");

    /* Get arguments */
    struct guppi_thread_args *args = (struct guppi_thread_args *)_args;

    /* Set cpu affinity */
    int rv = sched_setaffinity(0, sizeof(cpu_set_t), &args->cpuset);
    if (rv<0) { 
        guppi_error("guppi_pktsock_thread_codd", "Error setting cpu affinity.");
        perror("sched_setaffinity");
    }

    /* Set priority */
    // rv = setpriority(PRIO_PROCESS, 0, args->priority);
    rv=0;
    if (args->priority != 0)
    {
        struct sched_param priority_param;
        priority_param.sched_priority = args->priority;
        rv = pthread_setschedparam(pthread_self(), SCHED_FIFO, &priority_param);
    }    
    
    if (rv!=0) {
        guppi_error("guppi_pktsock_thread_codd", "Error setting priority level.");
        perror("set_priority");
    }
    printf("DBUG: codd_pktsockt thread pid=%d\n", gettid());

    /* Attach to status shared mem area */
    struct guppi_status st;
    rv = guppi_status_attach(&st);
    if (rv!=GUPPI_OK) {
        guppi_error("guppi_pktsock_thread_codd", 
                "Error attaching to status shared memory.");
        pthread_exit(NULL);
    }
    pthread_cleanup_push((void *)guppi_status_detach, &st);
    pthread_cleanup_push((void *)set_exit_status, &st);

    /* Init status, read info */
    guppi_status_lock_safe(&st);
    hputs(st.buf, STATUS_KEY, "init");
    hputs(st.buf, "DAQSTATE", "armed");
    guppi_status_unlock_safe(&st);

    /* Read in general parameters */
    struct guppi_params gp;
    struct psrfits pf;
    pf.sub.dat_freqs = NULL;
    pf.sub.dat_weights = NULL;
    pf.sub.dat_offsets = NULL;
    pf.sub.dat_scales = NULL;
    char status_buf[GUPPI_STATUS_SIZE];
    guppi_status_lock_safe(&st);
    memcpy(status_buf, st.buf, GUPPI_STATUS_SIZE);
    guppi_status_unlock_safe(&st);
    guppi_read_obs_params(status_buf, &gp, &pf);
    pthread_cleanup_push((void *)guppi_free_psrfits, &pf);

    /* Read network params */
    struct guppi_pktsock_params ps_params;
    guppi_read_pktsock_params(status_buf, &ps_params);

    /* Attach to databuf shared mem */
    struct guppi_databuf *db;
    db = guppi_databuf_attach(args->output_buffer); 
    if (db==NULL) {
        guppi_error("guppi_pktsock_thread_codd",
                "Error attaching to databuf shared memory.");
        pthread_exit(NULL);
    }
    pthread_cleanup_push((void *)guppi_databuf_detach, db);

    /* Set up packet socket */

    // Make frame_size be a divisor of block size so that frames will be
    // contiguous in mapped mempory.  block_size must also be a multiple of
    // page_size.  Easiest way is to oversize the frames to be 16384 bytes, which
    // is bigger than we need, but keeps things easy.
    ps_params.ps.frame_size = PKTSOCK_BYTES_PER_FRAME;
    // total number of frames
    ps_params.ps.nframes = PKTSOCK_NFRAMES;
    // number of blocks
    ps_params.ps.nblocks = PKTSOCK_NBLOCKS;

    int rv = hashpipe_pktsock_open(&ps_params.ps, ps_params.ifname, PACKET_RX_RING);
    if (rv!=HASHPIPE_OK) {
        guppi_error("guppi_pktsock_thread_codd", "Error opening pktsock.");
        pthread_exit(NULL);
    }
    pthread_cleanup_push((void (*)(void *))hashpipe_pktsock_close, &ps_params.ps);

    /* Time parameters */
    int stt_imjd=0, stt_smjd=0;
    double stt_offs=0.0;

    /* See which packet format to use */
    int use_parkes_packets=0, baseband_packets=1;
    int nchan=0, npol=0, acclen=0;
    nchan = pf.hdr.nchan;
    npol = pf.hdr.npol;
    if (strncmp(ps_params.packet_format, "PARKES", 6)==0) { use_parkes_packets=1; }
    if (use_parkes_packets) {
        printf("guppi_pktsock_thread_codd: Using Parkes UDP packet format.\n");
        acclen = gp.decimation_factor;
        if (acclen==0) { 
            guppi_error("guppi_pktsock_thread_codd", 
                    "ACC_LEN must be set to use Parkes format");
            pthread_exit(NULL);
        }
    }

    /* Figure out size of data in each packet, number of packets
     * per block, etc.  Changing packet size during an obs is not
     * recommended.
     */
    int block_size;
    size_t packet_data_size = guppi_udp_packet_datasize(ps_params.packet_size);
    if (use_parkes_packets) 
        packet_data_size = parkes_udp_packet_datasize(ps_params.packet_size);
    unsigned packets_per_block; 
    if (hgeti4(status_buf, "BLOCSIZE", &block_size)==0) {
            block_size = db->block_size;
            hputi4(status_buf, "BLOCSIZE", block_size);
    } else {
        if (block_size > db->block_size) {
            guppi_error("guppi_pktsock_thread_codd", "BLOCSIZE > databuf block_size");
            block_size = db->block_size;
            hputi4(status_buf, "BLOCSIZE", block_size);
        }
    }
    packets_per_block = block_size / packet_data_size;

    const uint64_t samples_per_second = 3UL*1000*1000*1000; // 3 GHz
    const uint64_t samples_per_spectrum = 1024;
    const uint64_t spectra_per_packet = 32;

    double dwell_seconds = 0.0;
    uint64_t dwell_blocks = 0;

    /* If we're in baseband mode, figure out how much to overlap
     * the data blocks.
     */
    int overlap_packets=0;
    if (baseband_packets) {
        if (hgeti4(status_buf, "OVERLAP", &overlap_packets)==0) {
            overlap_packets = 0; // Default to no overlap
        } else {
            // XXX This is only true for 8-bit, 2-pol data:
            int samples_per_packet = packet_data_size / nchan / (size_t)4;
            if (overlap_packets % samples_per_packet) {
                guppi_error("guppi_pktsock_thread_codd", 
                        "Overlap is not an integer number of packets");
                overlap_packets = (overlap_packets/samples_per_packet+1);
                hputi4(status_buf, "OVERLAP", 
                        overlap_packets*samples_per_packet);
            } else {
                overlap_packets = overlap_packets/samples_per_packet;
            }
        }
    }

    /* List of databuf blocks currently in use */
    unsigned i;
    const int nblock = 2;
    struct datablock_stats blocks[nblock];
    for (i=0; i<nblock; i++) 
        init_block(&blocks[i], db, packet_data_size, packets_per_block, 
                overlap_packets);

    /* Convenience names for first/last blocks in set */
    struct datablock_stats *fblock, *lblock;
    fblock = &blocks[0];
    lblock = &blocks[nblock-1];

    /* Misc counters, etc */
    char *curdata=NULL, *curheader=NULL;
    unsigned long long start_seq_num=0, stop_seq_num=0, seq_num, last_seq_num=2048, nextblock_seq_num=0;
    long long seq_num_diff;
    double drop_frac_avg=0.0;
    const double drop_lpf = 0.25;
    int netbuf_full = 0;
    char netbuf_status[128] = {};

    // Drop all packets to date
    unsigned char *p_frame;
    while((p_frame=hashpipe_pktsock_recv_frame_nonblock(&ps_params.ps))) {
        hashpipe_pktsock_release_frame(p_frame);
    }

    guppi_warn("guppi_pktsock_thread_codd", "Thread started!");

    /* Main loop */
    unsigned force_new_block=0, waiting=-1;
    signal(SIGINT,cc);
    while (run) {

        /* Wait for data */
        do {
            p_frame = hashpipe_pktsock_recv_udp_frame(
                &ps_params.ps, ps_params.port, 1000); // 1 second timeout

            /* Set "waiting" flag */
            if (!p_frame && run && waiting!=1) {
                guppi_status_lock_safe(&st);
                hputs(st.buf, STATUS_KEY, "waiting");
                guppi_status_unlock_safe(&st);
                waiting=1;
            }
        } while (!p_frame && run);

        if(!run) {
            // We're outta here!
            if(p_frame) {
                hashpipe_pktsock_release_frame(p_frame);
            }
            break;
        }

        // If we somehow get here with a NULL p_frame, just redo the loop
        if(!p_frame) {
            continue;
        }

        /* Check packet size */
        if(ps_params.packet_size == 0) {
            ps_params.packet_size = PKT_UDP_SIZE(p_frame) - 8;
        } else if(ps_params.packet_size != PKT_UDP_SIZE(p_frame) - 8) {
            /* Unexpected packet size, ignore? */
            nbogus_total++;
            if(nbogus_total % 1000000 == 0) {
                guppi_status_lock_safe(&st);
                hputi4(st.buf, "NBOGUS", nbogus_total);
                hputi4(st.buf, "PKTSIZE", PKT_UDP_SIZE(p_frame)-8);
                guppi_status_unlock_safe(&st);
            }
            // Release frame!
            hashpipe_pktsock_release_frame(p_frame);
            continue; 
        }

        /* Update status if needed */
        if (waiting!=0) {
            guppi_status_lock_safe(&st);
            hputs(st.buf, STATUS_KEY, "receiving");
            guppi_status_unlock_safe(&st);
            waiting=0;
        }

        /* Convert packet format if needed */
        if (use_parkes_packets) {
            parkes_to_guppi_from_payload(
                (char *)PKT_UDP_DATA(p_frame), acclen, npol, nchan);
        }

        /* Check seq num diff */
        seq_num = guppi_pktsock_seq_num(p_frame);
        seq_num_diff = seq_num - last_seq_num;
        if (seq_num_diff<=0) { 
            if (seq_num_diff<-1024) {
              time_t now = time(NULL);
              char * ts = ctime(&now);
              ts[24] = '\0';
              guppi_warn(ts, "guppi_pktsock_thread_codd forcing new frame. (reset/resync?)");
              force_new_block=1;
            } else if (seq_num_diff==0) {
                char msg[256];
                sprintf(msg, "Received duplicate packet (seq_num=%lld)", 
                        seq_num);
                guppi_warn("guppi_pktsock_thread_codd", msg);
            }
            else {
              // Release frame!
              hashpipe_pktsock_release_frame(p_frame);
              /* No going backwards */
              continue;
            }
        } else { 
            force_new_block=0; 
            npacket_total += seq_num_diff;
            ndropped_total += seq_num_diff - 1;
        }
        last_seq_num = seq_num;

        /* Determine if we go to next block */
        if ((seq_num>=nextblock_seq_num) || force_new_block) {

            /* Update drop stats */
            if (fblock->npacket)  
                drop_frac_avg = (1.0-drop_lpf)*drop_frac_avg 
                    + drop_lpf * 
                    (double)fblock->ndropped / 
                    (double)fblock->npacket;

            guppi_status_lock_safe(&st);
            hputi4(st.buf, "PKTIDX", fblock->packet_idx);                          
            hputr8(st.buf, "DROPAVG", drop_frac_avg);
            hputr8(st.buf, "DROPTOT", 
                    npacket_total ? 
                    (double)ndropped_total/(double)npacket_total 
                    : 0.0);
            hputr8(st.buf, "DROPBLK", 
                    fblock->npacket ? 
                    (double)fblock->ndropped/(double)fblock->npacket
                    : 0.0);

            // Calculate stop packet index using current value of SCANLEN
            hgetr8(st.buf, "SCANLEN", &dwell_seconds);

            // Dwell blocks is equal to:
            //
            //           dwell_seconds * samples/second
            //     -------------------------------------------
            //     samples/spectrum * spectra/pkt * pkts/block
            //
            // To get an integer number of blocks, simply truncate
            dwell_blocks = trunc(dwell_seconds * samples_per_second /
                    (samples_per_spectrum * spectra_per_packet * packets_per_block));

            stop_seq_num = start_seq_num + packets_per_block * dwell_blocks;
            hputi8(st.buf, "PKTSTOP", stop_seq_num);

            guppi_status_unlock_safe(&st);

            /* Finalize first block, and push it off the list.
             * Then grab next available block.
             */
            if (fblock->block_idx>=0) finalize_block(fblock);
            block_stack_push(blocks, nblock);
            increment_block(lblock, seq_num);
            curdata = guppi_databuf_data(db, lblock->block_idx);
            curheader = guppi_databuf_header(db, lblock->block_idx);
            nextblock_seq_num = lblock->packet_idx 
                + packets_per_block - overlap_packets;

            /* If new obs started, reset total counters, get start
             * time.  Start time is rounded to nearest integer
             * second, with warning if we're off that by more
             * than 100ms.  Any current blocks on the stack
             * are also finalized/reset */
            if (force_new_block) {

                /* Reset stats */
                npacket_total=0;
                ndropped_total=0;
                nbogus_total=0;

                /* Get obs start time */
                get_current_mjd(&stt_imjd, &stt_smjd, &stt_offs);
                if (stt_offs>0.5) { stt_smjd+=1; stt_offs-=1.0; }
                if (fabs(stt_offs)>0.1) { 
                    char msg[256];
                    sprintf(msg, 
                            "Second fraction = %3.1f ms > +/-100 ms",
                            stt_offs*1e3);
                    guppi_warn("guppi_pktsock_thread_codd", msg);
                }
                stt_offs = 0.0;

                /* Warn if 1st packet number is not zero */
                if (seq_num!=0) {
                    char msg[256];
                    sprintf(msg, "First packet number is not 0 (seq_num=%lld)",
                            seq_num);
                    guppi_warn("guppi_pktsock_thread_codd", msg);
                }

                /* Flush any current buffers */
                for (i=0; i<nblock-1; i++) {
                    if (blocks[i].block_idx>=0) 
                        finalize_block(&blocks[i]);
                    reset_block(&blocks[i]);
                }

            }

            /* Read/update current status shared mem */
            guppi_status_lock_safe(&st);
            if(seq_num >= stop_seq_num) {
                hputs(st.buf, "DAQSTATE", "armed");
                hputi4(st.buf, "STTVALID", 0);
            } else if (stt_imjd!=0) {
#if 0
                hputi4(st.buf, "STT_IMJD", stt_imjd);
                hputi4(st.buf, "STT_SMJD", stt_smjd);
#endif
                hputr8(st.buf, "STT_OFFS", stt_offs);
                int sttvalid;
                hgeti4(st.buf, "STTVALID", &sttvalid);
                if(!sttvalid) {
                    time_t now = time(NULL);
                    char * ts = ctime(&now);
                    ts[24] = '\0';
                    guppi_warn(ts, "guppi_pktsock_thread_codd STTVALID is 0, forcing STTVALID=1");
                }
                hputs(st.buf, "DAQSTATE", "record");
                hputi4(st.buf, "STTVALID", 1);
            }
#if 0
            else {
                // Put a non-accurate start time to avoid polyco 
                // errors.
                get_current_mjd(&stt_imjd, &stt_smjd, &stt_offs);
                hputi4(st.buf, "STT_IMJD", stt_imjd);
                hputi4(st.buf, "STT_SMJD", stt_smjd);
                hputi4(st.buf, "STTVALID", 0);
                // Reset to zero
                stt_imjd = 0;
                stt_smjd = 0;
            }
#endif
            memcpy(status_buf, st.buf, GUPPI_STATUS_SIZE);
            guppi_status_unlock_safe(&st);

            /* block size possibly changed on new obs */
            /* TODO: what about overlap...
             * Also, should this even be allowed ?
             */
            if (force_new_block) {
                if (hgeti4(status_buf, "BLOCSIZE", &block_size)==0) {
                        block_size = db->block_size;
                } else {
                    if (block_size > db->block_size) {
                        guppi_error("guppi_pktsock_thread_codd", 
                                "BLOCSIZE > databuf block_size");
                        block_size = db->block_size;
                    }
                }
                packets_per_block = block_size / packet_data_size;
            }
            hputi4(status_buf, "BLOCSIZE", block_size);

            /* Wait for new block to be free, then clear it
             * if necessary and fill its header with new values.
             */
            netbuf_full = guppi_databuf_total_status(db);
            sprintf(netbuf_status, "%d/%d", netbuf_full, db->n_block);
            guppi_status_lock_safe(&st);
            hputs(st.buf, STATUS_KEY, "waitfree");
            hputs(st.buf, "NETBUFST", netbuf_status);
            guppi_status_unlock_safe(&st);
            while ((rv=guppi_databuf_wait_free(db, lblock->block_idx)) 
                    != GUPPI_OK) {
                if (rv==GUPPI_TIMEOUT) {
                    waiting=1;
                    netbuf_full = guppi_databuf_total_status(db);
                    sprintf(netbuf_status, "%d/%d", netbuf_full, db->n_block);
                    guppi_status_lock_safe(&st);
                    hputs(st.buf, STATUS_KEY, "blocked");
                    hputs(st.buf, "NETBUFST", netbuf_status);
                    guppi_status_unlock_safe(&st);
                    continue;
                } else {
                    guppi_error("guppi_pktsock_thread_codd", 
                            "error waiting for free databuf");
                    run=0;
                    pthread_exit(NULL);
                    break;
                }
            }
            guppi_status_lock_safe(&st);
            hputs(st.buf, STATUS_KEY, "receiving");
            guppi_status_unlock_safe(&st);

            memcpy(curheader, status_buf, GUPPI_STATUS_SIZE);
            //if (baseband_packets) { memset(curdata, 0, block_size); }
            if (1) { memset(curdata, 0, block_size); }

        }

        /* Copy packet into any blocks where it belongs.
         * The "write packets" functions also update drop stats 
         * for blocks, etc.
         */
        for (i=0; i<nblock; i++) {
            if ((blocks[i].block_idx>=0) 
                    && (block_packet_check(&blocks[i],seq_num)==0)) {
                if (baseband_packets) 
                    write_baseband_packet_to_block_from_pktsock_frame(
                        &blocks[i], p_frame, nchan);
                else
                    write_search_packet_to_block_from_pktsock_frame(
                        &blocks[i], p_frame);
            }
        }

        // Release frame back to ring buffer
        hashpipe_pktsock_release_frame(p_frame);

        /* Will exit if thread has been cancelled */
        pthread_testcancel();
    }

    guppi_warn("guppi_pktsock_thread_codd", "Thread exiting!");
    pthread_exit(NULL);

    /* Have to close all push's */
    pthread_cleanup_pop(0); /* Closes push(hashpipe_pktsock_close) */
    pthread_cleanup_pop(0); /* Closes set_exit_status */
    pthread_cleanup_pop(0); /* Closes guppi_free_psrfits */
    pthread_cleanup_pop(0); /* Closes guppi_status_detach */
    pthread_cleanup_pop(0); /* Closes guppi_databuf_detach */
}

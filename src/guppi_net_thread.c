/* guppi_net_thread.c
 *
 * Routine to read packets from network and put them
 * into shared memory blocks.
 */

#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <sched.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>

#include "fitshead.h"
#include "write_psrfits.h"
#include "guppi_error.h"
#include "guppi_status.h"
#include "guppi_databuf.h"
#include "guppi_udp.h"
#include "guppi_time.h"

#define STATUS_KEY "NETSTAT"  /* Define before guppi_threads.h */
#include "guppi_threads.h"

/* Write info from param struct into fits-style buffer */
extern void guppi_write_params(char *buf, 
                               struct guppi_params *g,
                               struct psrfits *p);

/* Parse info from buffer into param struct */
extern void guppi_read_params(char *buf, 
                              struct guppi_params *g,
                              struct psrfits *p);

/* This thread is passed a single arg, pointer
 * to the guppi_udp_params struct.  This thread should 
 * be cancelled and restarted if any hardware params
 * change, as this potentially affects packet size, etc.
 */
void *guppi_net_thread(void *_up) {

    /* Set cpu affinity */
    cpu_set_t cpuset, cpuset_orig;
    sched_getaffinity(0, sizeof(cpu_set_t), &cpuset_orig);
    CPU_ZERO(&cpuset);
    CPU_SET(3, &cpuset);
    int rv = sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
    if (rv<0) { 
        guppi_error("guppi_net_thread", "Error setting cpu affinity.");
        perror("sched_setaffinity");
    }

    /* Set priority */
    rv = setpriority(PRIO_PROCESS, 0, 0);
    if (rv<0) {
        guppi_error("guppi_net_thread", "Error setting priority level.");
        perror("set_priority");
    }

    /* Get UDP param struct */
    struct guppi_udp_params *up =
        (struct guppi_udp_params *)_up;

    /* Attach to status shared mem area */
    struct guppi_status st;
    rv = guppi_status_attach(&st);
    if (rv!=GUPPI_OK) {
        guppi_error("guppi_net_thread", 
                "Error attaching to status shared memory.");
        pthread_exit(NULL);
    }
    pthread_cleanup_push((void *)set_exit_status, &st);

    /* Init status */
    guppi_status_lock_safe(&st);
    hputs(st.buf, STATUS_KEY, "init");
    guppi_status_unlock_safe(&st);

    /* Read in general parameters */
    char status_buf[GUPPI_STATUS_SIZE];
    struct guppi_params gp;
    struct psrfits pf;
    guppi_status_lock_safe(&st);
    guppi_read_params(st.buf, &gp, &pf);
    memcpy(status_buf, st.buf, GUPPI_STATUS_SIZE);
    guppi_status_unlock_safe(&st);

    /* Attach to databuf shared mem */
    struct guppi_databuf *db;
    db = guppi_databuf_attach(1); // TODO don't hardcode this 1
    if (db==NULL) {
        guppi_error("guppi_net_thread",
                "Error attaching to databuf shared memory.");
        pthread_exit(NULL);
    }

    /* Set up UDP socket */
    rv = guppi_udp_init(up);
    if (rv!=GUPPI_OK) {
        guppi_error("guppi_net_thread",
                "Error opening UDP socket.");
        pthread_exit(NULL);
    }
    pthread_cleanup_push((void *)guppi_udp_close, up);

    /* Time parameters */
    int stt_imjd=0, stt_smjd=0;
    double stt_offs=0.0;

    /* Figure out size of data in each packet, number of packets
     * per block, etc.
     * TODO : Figure out how/if to deal with packet size changing.
     */
    struct guppi_udp_packet p;
    size_t packet_data_size = guppi_udp_packet_datasize(up->packet_size); 
    unsigned packets_per_block = db->block_size / packet_data_size;

    /* Counters */
    unsigned long long npacket_total=0, npacket_block=0;
    unsigned long long ndropped_total=0, ndropped_block=0;
    unsigned long long nbogus_total=0, nbogus_block=0;
    unsigned long long curblock_seq_num=0, nextblock_seq_num=0;
    unsigned long long seq_num, last_seq_num=2048;
    int curblock=-1;
    char *curheader=NULL, *curdata=NULL;
    unsigned block_packet_idx=0, last_block_packet_idx=0;
    double drop_frac_avg=0.0;
    const double drop_lpf = 0.25;

    /* Main loop */
    unsigned i, force_new_block=0, waiting=-1;
    char *dataptr;
    long long seq_num_diff;
    signal(SIGINT,cc);
    while (run) {

        /* Wait for data */
        rv = guppi_udp_wait(up);
        if (rv!=GUPPI_OK) {
            if (rv==GUPPI_TIMEOUT) { 
                /* Set "waiting" flag */
                if (waiting!=1) {
                    guppi_status_lock_safe(&st);
                    hputs(st.buf, STATUS_KEY, "waiting");
                    guppi_status_unlock_safe(&st);
                    waiting=1;
                }
                continue; 
            } else {
                guppi_error("guppi_net_thread", 
                        "guppi_udp_wait returned error");
                perror("guppi_udp_wait");
                pthread_exit(NULL);
            }
        }

        /* Read packet */
        rv = guppi_udp_recv(up, &p);
        if (rv!=GUPPI_OK) {
            if (rv==GUPPI_ERR_PACKET) {
                /* Unexpected packet size, ignore? */
                nbogus_total++;
                nbogus_block++;
                continue; 
            } else {
                guppi_error("guppi_net_thread", 
                        "guppi_udp_recv returned error");
                perror("guppi_udp_recv");
                pthread_exit(NULL);
            }
        }

        /* Update status if needed */
        if (waiting!=0) {
            guppi_status_lock_safe(&st);
            hputs(st.buf, STATUS_KEY, "receiving");
            guppi_status_unlock_safe(&st);
            waiting=0;
        }

        /* Check seq num diff */
        seq_num = guppi_udp_packet_seq_num(&p);
        seq_num_diff = seq_num - last_seq_num;
        if (seq_num_diff<0) { 
            if (seq_num_diff<-1024) { force_new_block=1; }
            else  { continue; } /* No going backwards */
        } else { force_new_block=0; }

        /* Determine if we go to next block */
        if ((seq_num>=nextblock_seq_num) || force_new_block) {

            if (curblock>=0) { 
                /* Close out current block */
                hputi4(curheader, "PKTIDX", curblock_seq_num);
                hputi4(curheader, "PKTSIZE", packet_data_size);
                hputi4(curheader, "NPKT", npacket_block);
                hputi4(curheader, "NDROP", ndropped_block);
                guppi_write_params(curheader, &gp, &pf);
                guppi_databuf_set_filled(db, curblock);
            }

            if (npacket_block) { 
                drop_frac_avg = (1.0-drop_lpf)*drop_frac_avg 
                    + drop_lpf*(double)ndropped_block/(double)npacket_block;
            }

            /* Put drop stats in general status area */
            hputr8(st.buf, "DROPAVG", drop_frac_avg);
            hputr8(st.buf, "DROPTOT", 
                    npacket_total ? 
                    (double)ndropped_total/(double)npacket_total 
                    : 0.0);
            hputr8(st.buf, "DROPBLK", 
                    npacket_block ? 
                    (double)ndropped_block/(double)npacket_block 
                    : 0.0);

            /* Reset block counters */
            npacket_block=0;
            ndropped_block=0;
            nbogus_block=0;

            /* If new obs started, reset total counters, get start
             * time. */
            if (force_new_block) {
                npacket_total=0;
                ndropped_total=0;
                nbogus_total=0;
                get_current_mjd(&stt_imjd, &stt_smjd, &stt_offs);
            }

            /* Read/update current status shared mem */
            guppi_status_lock_safe(&st);
            guppi_read_params(st.buf, &gp, &pf);
            hputi4(st.buf, "STT_IMJD", stt_imjd);
            hputi4(st.buf, "STT_SMJD", stt_smjd);
            hputr8(st.buf, "STT_OFFS", stt_offs);
            memcpy(status_buf, st.buf, GUPPI_STATUS_SIZE);
            guppi_status_unlock_safe(&st);

            /* Advance to next block when free, update its header */
            curblock = (curblock + 1) % db->n_block;
            curheader = guppi_databuf_header(db, curblock);
            curdata = guppi_databuf_data(db, curblock);
            last_block_packet_idx = 0;
            curblock_seq_num = seq_num - (seq_num % packets_per_block);
            nextblock_seq_num = curblock_seq_num + packets_per_block;
            guppi_databuf_wait_free(db, curblock);
            memcpy(curheader, status_buf, GUPPI_STATUS_SIZE);
        }

        /* Skip dropped blocks, put packet in right spot */
        block_packet_idx = seq_num - curblock_seq_num;
        dataptr = curdata + last_block_packet_idx*packet_data_size;
        for (i=last_block_packet_idx; i<block_packet_idx; i++) {
            memset(dataptr, 0, packet_data_size);
            dataptr += packet_data_size;
            ndropped_block++;
            ndropped_total++;
            npacket_total++;
            npacket_block++;
        }
        memcpy(dataptr, guppi_udp_packet_data(&p), packet_data_size);
        npacket_total++;
        npacket_block++;
        last_block_packet_idx = block_packet_idx + 1;
        last_seq_num = seq_num;

        /* Will exit if thread has been cancelled */
        pthread_testcancel();
    }

    pthread_exit(NULL);

    /* Have to close all push's */
    pthread_cleanup_pop(0); /* Closes push(guppi_udp_close) */
    pthread_cleanup_pop(0); /* Closes set_exit_status */
}

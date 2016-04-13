/* guppi_rawdisk_thread.c
 *
 * Write databuf blocks out to disk.
 */

#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <sched.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "fitshead.h"
#include "psrfits.h"
#include "quantization.h"
#include "guppi_error.h"
#include "guppi_params.h"
#include "guppi_status.h"
#include "guppi_databuf.h"

#define STATUS_KEY "DISKSTAT"
#include "guppi_threads.h"

// 80 character string for the BACKEND header record.
static const char BACKEND_RECORD[] =
// 0000000000111111111122222222223333333333
// 0123456789012345678901234567890123456789
  "BACKEND = 'GUPPI   '                    " \
  "                                        ";

// Read a status buffer all of the key observation paramters
extern void guppi_read_obs_params(char *buf, 
                                  struct guppi_params *g, 
                                  struct psrfits *p);

/* Parse info from buffer into param struct */
extern void guppi_read_subint_params(char *buf, 
                                     struct guppi_params *g,
                                     struct psrfits *p);

ssize_t write_all(int fd, const void *buf, size_t bytes_to_write)
{
  size_t bytes_remaining = bytes_to_write;
  ssize_t bytes_written = 0;
  while(bytes_remaining != 0) {
    bytes_written = write(fd, buf, bytes_remaining);
    if(bytes_written == -1) {
      // Error!
      return -1;
    }
    bytes_remaining -= bytes_written;
    buf += bytes_written;
  }
  // All done!
  return bytes_to_write;
}

int safe_close(int *pfd) {
    if (pfd==NULL) return 0;
    fsync(*pfd);
    return close(*pfd);
}

void guppi_rawdisk_thread(void *_args) {
    /* Get args */
    struct guppi_thread_args *args = (struct guppi_thread_args *)_args;
    
    /* Set cpu affinity */
    int rv = sched_setaffinity(0, sizeof(cpu_set_t), &args->cpuset);
    if (rv<0) { 
        guppi_error("guppi_rawdisk_thread", "Error setting cpu affinity.");
        perror("sched_setaffinity");
    }


    /* Set priority */
    // rv = setpriority(PRIO_PROCESS, 0, 0);
    rv=0;
    if (args->priority != 0)
    {
        struct sched_param priority_param;
        priority_param.sched_priority = args->priority;
        rv = pthread_setschedparam(pthread_self(), SCHED_FIFO, &priority_param);
    }    
    
    if (rv!=0) {
        guppi_error("guppi_rawdisk_thread", "Error setting priority level.");
        perror("set_priority");
    }

    /* Attach to status shared mem area */
    struct guppi_status st;
    rv = guppi_status_attach(&st);
    if (rv!=GUPPI_OK) {
        guppi_error("guppi_rawdisk_thread", 
                "Error attaching to status shared memory.");
        pthread_exit(NULL);
    }
    pthread_cleanup_push((void *)guppi_status_detach, &st);
    pthread_cleanup_push((void *)set_exit_status, &st);

    /* Init status */
    guppi_status_lock_safe(&st);
    hputs(st.buf, STATUS_KEY, "init");
    guppi_status_unlock_safe(&st);

    /* Read in general parameters */
    struct guppi_params gp;
    struct psrfits pf;
    pf.sub.dat_freqs = NULL;
    pf.sub.dat_weights = NULL;
    pf.sub.dat_offsets = NULL;
    pf.sub.dat_scales = NULL;
    pthread_cleanup_push((void *)guppi_free_psrfits, &pf);

    /* Attach to databuf shared mem */
    struct guppi_databuf *db;
    db = guppi_databuf_attach(args->input_buffer);
    if (db==NULL) {
        guppi_error("guppi_rawdisk_thread",
                "Error attaching to databuf shared memory.");
        pthread_exit(NULL);
    }
    pthread_cleanup_push((void *)guppi_databuf_detach, db);

    /* Init output file */
    static int fdraw = 0;
    pthread_cleanup_push((void *)safe_close, &fdraw);

    /* Pointers for quantization params */
    //double *mean = NULL;
    //double *std = NULL;
    double mean[2*256];
    double std[2*256];

    /* Loop */
    int packetidx=0, npacket=0, ndrop=0, packetsize=0, blocksize=0, len=0;
    int orig_blocksize=0;
    int curblock=0;
    int block_count=0, blocks_per_file=128, filenum=0;
    int got_packet_0=0, first=1;
    int first_block = 0;
    int requantize = 0;
    char *ptr, *hend;
    int open_flags = 0;
    int directio = 0;
    signal(SIGINT,cc);
    while (run) {

        /* Note waiting status */
        guppi_status_lock_safe(&st);
        hputs(st.buf, STATUS_KEY, "waiting");
        guppi_status_unlock_safe(&st);

        /* Wait for buf to have data */
        rv = guppi_databuf_wait_filled(db, curblock);
        if (rv!=0) continue;

        /* Read param struct for this block */
        ptr = guppi_databuf_header(db, curblock);
        if (first) {
            guppi_read_obs_params(ptr, &gp, &pf);
            first = 0;
        } else {
            guppi_read_subint_params(ptr, &gp, &pf);
        }

        /* Parse packet size, npacket from header */
        hgeti4(ptr, "PKTIDX", &packetidx);
        hgeti4(ptr, "PKTSIZE", &packetsize);
        hgeti4(ptr, "NPKT", &npacket);
        hgeti4(ptr, "NDROP", &ndrop);

        /* Check for re-quantization flag */
        int nbits_req = 0;
        if (hgeti4(ptr, "NBITSREQ", &nbits_req) == 0) {
            /* Param not present, don't requantize */
            requantize = 0;
        } else {
            /* Param is present */
            if (nbits_req==8)
                requantize = 0;
            else if (nbits_req==2) 
                requantize = 1;
            else
                /* Invalid selection for requested nbits 
                 * .. die or ignore?
                 */
                requantize = 0;
        }

        /* Set up data ptr for quant routines */
        pf.sub.data = (unsigned char *)guppi_databuf_data(db, curblock);

        /* Wait for packet 0 before starting write */
        if (got_packet_0==0 && packetidx==0 && gp.stt_valid==1) {
            got_packet_0 = 1;
            guppi_read_obs_params(ptr, &gp, &pf);
            orig_blocksize = pf.sub.bytes_per_subint;
            directio = guppi_read_directio_mode(ptr);
            char fname[256];
            sprintf(fname, "%s.%4.4d.raw", pf.basefilename, filenum);
            fprintf(stderr, "Opening raw file '%s' (directio=%d)\n", fname, directio);
            // Create the output directory if needed
            char datadir[1024];
            strncpy(datadir, pf.basefilename, 1023);
            char *last_slash = strrchr(datadir, '/');
            if (last_slash!=NULL && last_slash!=datadir) {
                *last_slash = '\0';
                printf("Using directory '%s' for output.\n", datadir);
                char cmd[1024];
                sprintf(cmd, "mkdir -m 1777 -p %s", datadir);
                system(cmd);
            }
            // TODO: check for file exist.
            open_flags = O_CREAT|O_RDWR|O_SYNC;
            if(directio) {
              open_flags |= O_DIRECT;
            }
            fdraw = open(fname, open_flags, 0644);
            if (fdraw==-1) {
                guppi_error("guppi_rawdisk_thread", "Error opening file.");
                pthread_exit(NULL);
            }
            /* Determine scaling factors for quantization if appropriate */
            if (requantize) {
#if 0 
                mean = (double *)realloc(mean, 
                        pf.hdr.rcvr_polns * pf.hdr.nchan * sizeof(double));
                std  = (double *)realloc(std,  
                        pf.hdr.rcvr_polns * pf.hdr.nchan * sizeof(double));
                fprintf(stderr, "Alloced 2-bit arrays\n"); fflush(stderr);
                if (mean==NULL) {
                    fprintf(stderr, "mean is null\n"); fflush(stderr);
                }
                if (std==NULL) {
                    fprintf(stderr, "std is null\n"); fflush(stderr);
                }
#endif
                fprintf(stderr, "curblock=%d\n", curblock); fflush(stderr);
                compute_stat(&pf, mean, std);
                fprintf(stderr, "Computed 2-bit stats\n"); fflush(stderr);
                first_block = 1;
            }
        }
        
        /* See if we need to open next file */
        if (block_count >= blocks_per_file) {
            close(fdraw);
            filenum++;
            char fname[256];
            sprintf(fname, "%s.%4.4d.raw", pf.basefilename, filenum);
            directio = guppi_read_directio_mode(ptr);
            open_flags = O_CREAT|O_RDWR|O_SYNC;
            if(directio) {
              open_flags |= O_DIRECT;
            }
            fprintf(stderr, "Opening raw file '%s' (directio=%d)\n", fname, directio);
            fdraw = open(fname, open_flags, 0644);
            if (fdraw==-1) {
                guppi_error("guppi_rawdisk_thread", "Error opening file.");
                pthread_exit(NULL);
            }
            block_count=0;
        }

        /* See how full databuf is */
        //total_status = guppi_databuf_total_status(db);

        /* Requantize from 8 bits to 2 bits if necessary.
         * See raw_quant.c for more usage examples.
         */
        if (requantize && got_packet_0 && !first_block) {
        //if (requantize && got_packet_0) {
            pf.sub.bytes_per_subint = orig_blocksize;
            /* Does the quantization in-place */
            quantize_2bit(&pf, mean, std);
            /* Update some parameters for output */
            hputi4(ptr, "BLOCSIZE", pf.sub.bytes_per_subint);
            hputi4(ptr, "NBITS", pf.hdr.nbits);
        }

        /* Get full data block size */
        hgeti4(ptr, "BLOCSIZE", &blocksize);

        /* If we got packet 0, write data to disk */
        if (got_packet_0) { 

            /* Note waiting status */
            guppi_status_lock_safe(&st);
            hputs(st.buf, STATUS_KEY, "writing");
            guppi_status_unlock_safe(&st);

            /* Write header to file */
            hend = ksearch(ptr, "END");
            len = (hend-ptr)+80;

            // If BACKEND record is not present, insert it as first record.
            // TODO: Verify that we have room to insert the record.
            if(!ksearch(ptr, "BACKEND")) {
                // Move exsiting records to make room for new first record
                memmove(ptr+80, ptr, len);
                // Copy in BACKEND_RECORD string
                strncpy(ptr, BACKEND_RECORD, 80);
                // Increase len by 80 to account for the added record
                len += 80;
            }

            // Adjust length for any padding required for DirectIO
            if(directio) {
                // Round up to next multiple of 512
                len = (len+511) & ~511;
            }

            /* Write header (and padding, if any) */
            rv = write_all(fdraw, ptr, len);
            if (rv != len) {
                char msg[100];
                perror("guppi_rawdisk_thread write_all header");
                sprintf(msg, "Error writing data (ptr=%p, len=%d, rv=%d)", ptr, len, rv);
                guppi_error("guppi_rawdisk_thread", msg);
                        //"Error writing data.");
            }

            /* Write data */
            ptr = guppi_databuf_data(db, curblock);
            len = blocksize;
            if(directio) {
                // Round up to next multiple of 512
                len = (len+511) & ~511;
            }
            rv = write_all(fdraw, ptr, (size_t)len);
            if (rv != len) {
                char msg[100];
                perror("guppi_rawdisk_thread write_all block");
                sprintf(msg, "Error writing data (ptr=%p, len=%d, rv=%d)", ptr, len, rv);
                guppi_error("guppi_rawdisk_thread", msg);
                        //"Error writing data.");
            }

            /* Increment counter */
            block_count++;
            first_block = 0;

            /* flush output */
            fsync(fdraw);
        }

        /* Mark as free */
        guppi_databuf_set_free(db, curblock);

        /* Go to next block */
        curblock = (curblock + 1) % db->n_block;

        /* Check for cancel */
        pthread_testcancel();

    }

    pthread_exit(NULL);

    pthread_cleanup_pop(0); /* Closes close */
    pthread_cleanup_pop(0); /* Closes guppi_databuf_detach */
    pthread_cleanup_pop(0); /* Closes guppi_free_psrfits */
    pthread_cleanup_pop(0); /* Closes set_exit_status */
    pthread_cleanup_pop(0); /* Closes guppi_status_detach */

}

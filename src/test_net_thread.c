/* test_net_thread.c
 *
 * Test run net thread.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <signal.h>
#include <poll.h>
#include <getopt.h>
#include <errno.h>

#include "write_psrfits.h"
#include "guppi_udp.h"
#include "guppi_error.h"
#include "guppi_status.h"
#include "guppi_databuf.h"
#include "guppi_params.h"

#include "guppi_thread_main.h"

/* Write info from param struct into fits-style buffer */
extern void guppi_write_params(char *buf, 
                               struct guppi_params *g,
                               struct psrfits *p);

void usage() {
    fprintf(stderr,
            "Usage: test_net_thread [options] sender_hostname\n"
            "Options:\n"
            "  -p n, --port=n    Port number\n"
            "  -h, --help        This message\n"
           );
}

/* Thread declarations */
void *guppi_net_thread(void *_up);
void *guppi_rawdisk_thread(void *args);

int main(int argc, char *argv[]) {

    struct guppi_udp_params p;

    static struct option long_opts[] = {
        {"help",   0, NULL, 'h'},
        {"port",   1, NULL, 'p'},
        {0,0,0,0}
    };
    int opt, opti;
    p.port = 5000;
    while ((opt=getopt_long(argc,argv,"hp:",long_opts,&opti))!=-1) {
        switch (opt) {
            case 'p':
                p.port = atoi(optarg);
                break;
            default:
            case 'h':
                usage();
                exit(0);
                break;
        }
    }

    /* Need sender hostname */
    if (optind==argc) {
        usage();
        exit(1);
    }

    /* Init udp params */
    strcpy(p.sender, argv[optind]);
    p.packet_size = 8208; /* Expected 8k + 8 byte seq num + 8 byte flags */

    /* Init shared mem */
    struct guppi_status stat;
    struct guppi_databuf *dbuf=NULL;
    int rv = guppi_status_attach(&stat);
    if (rv!=GUPPI_OK) {
        fprintf(stderr, "Error connecting to guppi_status\n");
        exit(1);
    }
    dbuf = guppi_databuf_attach(1);
    if (dbuf==NULL) {
        fprintf(stderr, "Error connecting to guppi_databuf\n");
        exit(1);
    }
    guppi_databuf_clear(dbuf);

    /* Fake parameters */
    struct guppi_params gp;
    struct psrfits pf;
    pf.hdr.MJD_epoch = 54535.567;
    pf.hdr.fctr = 2000.0;
    pf.hdr.BW = 800.0;
    pf.hdr.nchan = 4096;
    pf.hdr.nbits = 8;
    pf.hdr.npol = 1;
    pf.hdr.dt = 81.92e-6;
    guppi_status_lock(&stat);
    guppi_write_params(stat.buf, &gp, &pf);
    guppi_status_unlock(&stat);

    run=1;
    signal(SIGINT, cc);

    /* Launch net thread */
    pthread_t net_thread_id;
    rv = pthread_create(&net_thread_id, NULL, guppi_net_thread,
            (void *)&p);
    if (rv) { 
        fprintf(stderr, "Error creating net thread.\n");
        perror("pthread_create");
        exit(1);
    }

    /* Launch raw disk thread */
    pthread_t disk_thread_id;
    rv = pthread_create(&disk_thread_id, NULL, guppi_rawdisk_thread, NULL);
    if (rv) { 
        fprintf(stderr, "Error creating net thread.\n");
        perror("pthread_create");
        exit(1);
    }

    /* Wait for end */
    while (run) { sleep(1); }
    pthread_cancel(disk_thread_id);
    pthread_cancel(net_thread_id);
    pthread_kill(disk_thread_id,SIGINT);
    pthread_kill(net_thread_id,SIGINT);
    pthread_join(net_thread_id,NULL);
    printf("Joined net thread\n"); fflush(stdout);
    pthread_join(disk_thread_id,NULL);
    printf("Joined disk thread\n"); fflush(stdout);

    exit(0);
}

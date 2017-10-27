/* check_guppi_databuf.c
 *
 * Basic prog to test dstabuf shared mem routines.
 */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/sem.h>
#include <errno.h>

#include "fitshead.h"
#include "guppi_error.h"
#include "guppi_status.h"
#include "guppi_databuf.h"

void usage() { 
    fprintf(stderr, 
            "Usage: check_guppi_databuf [options]\n"
            "Options:\n"
            "  -h, --help\n"
            "  -q, --quiet\n"
            "  -c, --create\n"
            "  -d, --delete\n"            
            "  -S, --status\n"
            "  -i n, --id=n  (1)\n"
            "  -s n, --size=n (32M)\n"
            "  -n n, --nblock=n (24)\n"
            );
}

int main(int argc, char *argv[]) {

    /* Loop over cmd line to fill in params */
    static struct option long_opts[] = {
        {"help",   0, NULL, 'h'},
        {"quiet",  0, NULL, 'q'},
        {"create", 0, NULL, 'c'},
        {"delete", 0, NULL, 'd'},        
        {"status", 0, NULL, 'S'},
        {"id",     1, NULL, 'i'},
        {"size",   1, NULL, 's'},
        {"nblock", 1, NULL, 'n'},
        {0,0,0,0}
    };
    int opt,opti;
    int quiet=0;
    int create=0;
    int db_id=1;
    int blocksize = 32;
    int nblock = 24;
    int deletebuf=0;
    int dostatus=0;
    while ((opt=getopt_long(argc,argv,"hqcdSi:s:n:",long_opts,&opti))!=-1) {
        switch (opt) {
            case 'c':
                create=1;
                break;
            case 'q':
                quiet=1;
                break;
            case 'i':
                db_id = atoi(optarg);
                break;
            case 's':
                blocksize = atoi(optarg);
                break;
            case 'n':
                nblock = atoi(optarg);
                break;
            case 'd':
                deletebuf=1;
                create=0;
                break;                
            case 'S':
                deletebuf=0;
                dostatus=1;
                create=0;
                break;
            case 'h':
            default:
                usage();
                exit(0);
                break;
        }
    }

    /* Create mem if asked, otherwise attach */
    struct guppi_databuf *db=NULL;
    if (create) { 
        db = guppi_databuf_create(nblock, blocksize*1024*1024, db_id);
        if (db==NULL) {
            fprintf(stderr, "Error creating databuf %d (may already exist).\n",
                    db_id);
            exit(1);
        }
    } else {
        db = guppi_databuf_attach(db_id);
        if (db==NULL) { 
            fprintf(stderr, 
                    "Error attaching to databuf %d (may not exist).\n",
                    db_id);
            exit(1);
        }
        if (deletebuf)
        {
            /* attach worked so it exists. Now clear it and detach in
               preparation to delete.
            */
            int shmid, semid, rtnval;
            guppi_databuf_clear(db);          
            shmid = db->shmid;
            semid = db->semid;
            rtnval = 0;
            if (shmctl(shmid,IPC_RMID, 0) != 0)
            {
                perror("removal of buffer failed:");
                rtnval = -1;
            }
            printf("buffer deleted successfully\n");
            if (shmdt(db) != 0)
            {
                perror("shm detach failed:");
                // silently fail
            }
            if (semctl(semid,IPC_RMID, 0) != 0)
            {
                perror("removal of semaphores failed:");
                rtnval = -1;
            }
            printf("sems deleted successfully\n");            
            exit (rtnval);
        }
        else if(dostatus)
        {
            size_t n = db->n_block + 1;
            char * s = malloc(n);
            if(s)
            {
                guppi_databuf_str_status(db, s, n);
                printf("%d %s\n", db_id, s);
                free(s);
            }
        }
    }
    
    if (quiet)
    {
        /* skip the verbose stats */
        exit(0);
    }

    /* Print basic info */
    printf("databuf %d stats:\n", db_id);
    printf("  shmid=%d\n", db->shmid);
    printf("  semid=%d\n", db->semid);
    printf("  n_block=%d\n", db->n_block);
    printf("  struct_size=%zd\n", db->struct_size);
    printf("  block_size=%zd\n", db->block_size);
    printf("  header_size=%zd\n\n", db->header_size);

    /* loop over blocks */
    int i;
    char buf[81];
    char *hdr, *ptr, *hend;
    for (i=0; i<db->n_block; i++) {
        printf("block %d status=%d\n", i, 
                guppi_databuf_block_status(db, i));
        hdr = guppi_databuf_header(db, i);
        hend = ksearch(hdr, "END");
        if (hend==NULL) {
            printf("header not initialized\n");
        } else {
            hend += 80;
            printf("header:\n");
            for (ptr=hdr; ptr<hend; ptr+=80) {
                strncpy(buf, ptr, 80);
                buf[79]='\0';
                printf("%s\n", buf);
            }
        }

    }

    exit(0);
}

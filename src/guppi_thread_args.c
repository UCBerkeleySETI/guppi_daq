#define _GNU_SOURCE 1

#include "guppi_thread_args.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/prctl.h>


void guppi_thread_args_init(struct guppi_thread_args *a) {
    a->priority=0;
    a->finished=0;
    pthread_cond_init(&a->finished_c,NULL);
    pthread_mutex_init(&a->finished_m,NULL);
    /* By default, allow all cores currently allowed */
    sched_getaffinity(0, sizeof(cpu_set_t), &a->cpuset);
}

void guppi_thread_args_destroy(struct guppi_thread_args *a) {
    a->finished=1;
    pthread_cond_destroy(&a->finished_c);
    pthread_mutex_destroy(&a->finished_m);
}

void guppi_thread_set_finished(struct guppi_thread_args *a) {
    pthread_mutex_lock(&a->finished_m);
    a->finished=1;
    pthread_cond_broadcast(&a->finished_c);
    pthread_mutex_unlock(&a->finished_m);
}

int guppi_thread_finished(struct guppi_thread_args *a, 
        float timeout_sec) {
    struct timeval now;
    struct timespec twait;
    int rv;
    pthread_mutex_lock(&a->finished_m);
    gettimeofday(&now,NULL);
    twait.tv_sec = now.tv_sec + (int)timeout_sec;
    twait.tv_nsec = now.tv_usec * 1000 + 
        (int)(1e9*(timeout_sec-floor(timeout_sec)));
    if (a->finished==0) 
        rv = pthread_cond_timedwait(&a->finished_c, &a->finished_m, &twait);
    rv = a->finished;
    pthread_mutex_unlock(&a->finished_m);
    return(rv);
}

// A thread affinity/thread priority config file reader
void read_thread_configuration(struct KeywordValues *keywords)
{
    char *conf_root = getenv("GUPPI_DIR");
    char conf_file_name[256];
    char linebuf[128];
    int i;
    FILE *fin;
    char *saveptr, *keyword, *value, *endptr;
    
    if (conf_root == 0)
    {
        printf("GUPPI_DIR is not set to config directory\n");
        return;
    }
    snprintf(conf_file_name, sizeof(conf_file_name), "%s/guppi_threads.conf", conf_root);
    fin = fopen(conf_file_name, "r");
    if (!fin)
    {
        printf("Warning: thread configuration file %s not found\n", conf_file_name);
        return;
    }
    do
    {
        if (fgets(linebuf, sizeof(linebuf), fin) == 0)
        {
            break;
        }
        // not a comment, not a blank line, and an '=' is present
        if (strlen(linebuf) > 3 && linebuf[0] != '#' && strchr(linebuf, '='))
        {
           int found_key;
           keyword = strtok_r(linebuf, "=" , &saveptr);
           value = strtok_r(NULL, "# " , &saveptr);
           for (found_key=0, i=0; !found_key && keywords[i].name; ++i)
           {
               if (keyword && !strcasecmp(keywords[i].name, keyword))
               {
                   unsigned int val = strtol(value, &endptr, 0);
                   if (value == endptr)
                   {
                       printf("Error reading numeric value on line:%s\n", linebuf);
                       continue;
                   }
                   keywords[i].value = val;
                   found_key=1;
               }
           }
        }
    } while (!feof(fin));
    fclose(fin);
}

unsigned int get_config_key_value(char *keyword, struct KeywordValues *keywords)
{
    int i;
    for (i=0; keywords[i].name; ++i)
    {
        if (!strcasecmp(keywords[i].name, keyword))
        {
            return keywords[i].value;
        }
    }
    return 0;
}

void mask_to_cpuset(cpu_set_t *cpuset, unsigned int mask)
{
    int i;
    int core=0;
    // illegal condtion -- use default
    if (mask == 0)
        return;
        
    CPU_ZERO(cpuset);
    for (i=1; i<31; ++i, ++core)
    {
        if ((1<<i) & mask)
        {
            CPU_SET(core, cpuset);
        }
    }
}


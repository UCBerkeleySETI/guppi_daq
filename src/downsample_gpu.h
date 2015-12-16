#ifndef _DOWNSAMPLE_GPU_H
#define _DOWNSAMPLE_GPU_H
#include "dedisperse_gpu.h"
#ifdef __cplusplus
extern "C" {
#endif
void init_downsample(struct dedispersion_setup *s);
void downsample(struct dedispersion_setup *s, char *ds_out);
void transpose8(struct dedispersion_setup *s,int nbytes, char *ds_out);
void transfer_to_host(struct dedispersion_setup *s, int big_ds_bytes, char *ds_out);
#ifdef __cplusplus
}
#endif
#endif

/* downsample_gpu.cu
 * Detection/downsampling in GPU/CUDA
 * Paul Demorest, 2009/10
 */
#include <math.h>
#include "dedisperse_gpu.h"
#include "downsample_gpu.h"

/* Returns number of bytes needed for downsampling block */
size_t get_ds_bytes(const struct dedispersion_setup *s) {
    if (s->npol==1) 
        return sizeof(char) * s->npts_per_block / s->dsfac;
    else
        return sizeof(char4) * s->npts_per_block / s->dsfac;
}

/* Initialize the downsampling using values in dedispersion_setup
 * struct.  the s->dsfac field needs to be filled in.
 */
extern "C"
void init_downsample(struct dedispersion_setup *s) {

    // TODO: check that params satisfy any alignment requirements.

    // Allocate memory for DS results on GPU
    const size_t ds_bytes = get_ds_bytes(s) * s->nchan;
    cudaMalloc((void**)&s->dsbuf_gpu, ds_bytes);
    cudaMalloc((void**)&s->dsbuf_trans_gpu, ds_bytes);
    printf("Downsample memory = %.1f MB\n", ds_bytes / (1024.*1024.));

    // Check for errors
    cudaThreadSynchronize();
    printf("init_downsample cuda_err='%s'\n", 
            cudaGetErrorString(cudaGetLastError()));
}

extern "C"
void deinit_downsample(struct dedispersion_setup *s)
{
    if (s->dsbuf_gpu)
    {
        cudaFree(s->dsbuf_gpu);
        s->dsbuf_gpu=0;
    }
    if (s->dsbuf_trans_gpu)
    {
        cudaFree(s->dsbuf_trans_gpu);
        s->dsbuf_trans_gpu = 0;
    }
    printf("deinit_downsample cuda_err='%s'\n", 
            cudaGetErrorString(cudaGetLastError()));
}

/* "naive" version where each thread does one output sample at a time
 * If this isn't fast enough there are lots of optimizations that
 * could be done...
 */
__global__ void detect_downsample_4pol(const float2 *pol0, const float2 *pol1,
        const unsigned dsfac, const unsigned fftlen, const unsigned overlap,
        char4 *out) {

    // Dimensions
    const int tid = threadIdx.x;
    const int nt = blockDim.x;
    const int nvalid = fftlen - overlap;
    const int ifft = blockIdx.x;
    const int iblock = blockIdx.y;
    const int nsamp_per_block = nvalid / gridDim.y;
    const int nout_per_block = nsamp_per_block / dsfac;

    // Data pointers
    const float2 *ptr0 = pol0 + ifft*fftlen + overlap/2 
        + iblock*nsamp_per_block;
    const float2 *ptr1 = pol1 + ifft*fftlen + overlap/2 
        + iblock*nsamp_per_block;
    char4 *optr = out + ifft*nvalid/dsfac + iblock*nout_per_block;

    // Data scaling
    // This should be appropriate for input baseband data with
    // a RMS of ~20 counts.
    const float scale = (float)dsfac * 20.0;

    // Loop over data
    for (int iout=tid; iout<nout_per_block; iout+=nt) {
        float4 otmp= make_float4(0,0,0,0);
        for (int j=0; j<dsfac; j++) {
            float2 p0 = ptr0[iout*dsfac+j];
            float2 p1 = ptr1[iout*dsfac+j];
            otmp.x += p0.x*p0.x + p0.y*p0.y;
            otmp.y += p1.x*p1.x + p1.y*p1.y;
            otmp.z += p0.x*p1.x + p0.y*p1.y;
            otmp.w += p0.x*p1.y - p0.y*p1.x;
        }
        // The __float2int_rn function appears to not saturate, so if otmp.x/scale > 255, it will wrap
        // so large spikes in the data (giant pulses!) could be lost
        // This code saturates to 255 so that the value 255 in the PP or QQ terms indicates that some 
        // clipping has occured, in which case the cross terms should not be trusted.
        otmp.x = otmp.x/scale;
        if (otmp.x > 254) {
            otmp.x = 255;
        }
        otmp.y = otmp.y/scale;
        if (otmp.y > 254) {
            otmp.y = 255;
        }
        optr[iout].x = __float2int_rn(otmp.x);
        optr[iout].y = __float2int_rn(otmp.y);
    
        optr[iout].z = __float2int_rn(otmp.z/scale);
        optr[iout].w = __float2int_rn(otmp.w/scale);
    }

}

/* Same as above, except only compute total power */
__global__ void detect_downsample_1pol(const float2 *pol0, const float2 *pol1,
        const unsigned dsfac, const unsigned fftlen, const unsigned overlap,
        char *out) {

    // Dimensions
    const int tid = threadIdx.x;
    const int nt = blockDim.x;
    const int nvalid = fftlen - overlap;
    const int ifft = blockIdx.x;
    const int iblock = blockIdx.y;
    const int nsamp_per_block = nvalid / gridDim.y;
    const int nout_per_block = nsamp_per_block / dsfac;

    // Data pointers
    const float2 *ptr0 = pol0 + ifft*fftlen + overlap/2 
        + iblock*nsamp_per_block;
    const float2 *ptr1 = pol1 + ifft*fftlen + overlap/2 
        + iblock*nsamp_per_block;
    char *optr = out + ifft*nvalid/dsfac + iblock*nout_per_block;

    // Data scaling
    // This should be appropriate for input baseband data with
    // a RMS of ~20 counts in each poln (final 2.0 is for polns).
    const float scale = (float)dsfac * 20.0 * 2.0;

    // Loop over data
    for (int iout=tid; iout<nout_per_block; iout+=nt) {
        float otmp = 0.0;
        for (int j=0; j<dsfac; j++) {
            float2 p0 = ptr0[iout*dsfac+j];
            float2 p1 = ptr1[iout*dsfac+j];
            otmp += p0.x*p0.x + p0.y*p0.y + p1.x*p1.x + p1.y*p1.y;
        }
        optr[iout] = __float2int_rn(otmp/scale);
    }

}
__global__
void gpu_transpose8(char *in, char *out)
{
    int dx = blockIdx.x * blockDim.x + threadIdx.x;
    int dy = blockIdx.y * blockDim.y + threadIdx.y;
    int in_idx  = dy * gridDim.x * blockDim.x + dx;

    // exchange the coordinates
    int out_idx = dx * gridDim.y * blockDim.y + dy;

    out[out_idx] = in[in_idx];
}

/// Perform a classic non-square transpose of byte values
/// @param ds_out - if NULL data is left on GPU, otherwise the data is
/// transferred back to the host in the pointer ds_out.
extern "C"
void transpose8(struct dedispersion_setup *s, int big_ds_bytes, char *ds_out)
{
    cudaEvent_t start, stop;
    float transp_time, copy_time;
    
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    dim3 grid, block;
    block.x = 8;
    block.y = 8;
    grid.x = ((s->npts_per_block/s->dsfac) * s->npol)/block.x;
    grid.y = (s->nchan/block.y);
    
    cudaEventRecord(start, 0);
    // CHECK_CUDA("transpose8(A)");
    
    gpu_transpose8<<<grid, block>>>(s->dsbuf_gpu, s->dsbuf_trans_gpu);
    CHECK_CUDA("transpose8(b)");
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&transp_time, start, stop);
    s->time.downsample += transp_time;
    s->time.total2 += transp_time;
    
    
    /* Transfer data back to CPU if ds_out is not null*/
    if (ds_out != 0)
    {
        cudaEventRecord(start, 0);
        cudaMemcpy(ds_out, s->dsbuf_trans_gpu, big_ds_bytes, cudaMemcpyDeviceToHost);
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&copy_time, start, stop);
        s->time.transfer_to_host += copy_time;
        s->time.total2 += copy_time;
    }
}

/** Detect / downsample data.  Assumes dedispersion results
 * are already in the GPU, as described in the dedispersion_setup
 * struct.
 * @param ds_out - if NULL data is left on GPU, otherwise the data is
 * transferred back to the host in the pointer ds_out.
 */
extern "C"
void downsample(struct dedispersion_setup *s, char *ds_out) {

    /* Sizes */
    const size_t ds_bytes = get_ds_bytes(s);

    /* Benchmark */
#define NT 5
    cudaEvent_t t[NT];
    int it;
    for (it=0; it<NT; it++) cudaEventCreate(&t[it]);
    it=0;

    cudaEventRecord(t[it], 0); it++;
    cudaEventRecord(t[it], 0); it++;
    // CHECK_CUDA("downsample(A)");
    /* Clear out data buf */
    cudaMemset(s->dsbuf_gpu, 0, ds_bytes);

    /* Downsample data */
    // Setup Grid dimensions: nffts_per_block x 32 [thread blocks]
    // Each thread block has 64 x 1 threads
    // What if we increased the dsbuf size and added a zdim of size nchan?
    // 
    // CHECK_CUDA("downsample(a)");
    dim3 gd(s->nfft_per_block, 32, 1);
    if (s->npol==1) 
        detect_downsample_1pol<<<gd, 64>>>(s->databuf0_gpu, s->databuf1_gpu,
                s->dsfac, s->fft_len, s->overlap, (char *)s->dsbuf_gpu);
    else
        detect_downsample_4pol<<<gd, 64>>>(s->databuf0_gpu, s->databuf1_gpu,
                s->dsfac, s->fft_len, s->overlap, (char4 *)s->dsbuf_gpu);
    cudaEventRecord(t[it], 0); it++;
    
    CHECK_CUDA("downsample(b)");

    /* Transfer data back to CPU if ds_out is not null*/
    if (ds_out != 0)
    {
        cudaMemcpy(ds_out, s->dsbuf_gpu, ds_bytes, cudaMemcpyDeviceToHost);
    }
    // CHECK_CUDA("downsample(c)");
    cudaEventRecord(t[it], 0); it++;

    /* Final timer */
    cudaEventRecord(t[it], 0);
    cudaEventSynchronize(t[it]);
    cudaThreadSynchronize();

    /* Add up timers */
    float ttmp;
    it=1;

    cudaEventElapsedTime(&ttmp, t[it], t[it+1]);
    s->time.downsample += ttmp;
    s->time.total2 += ttmp;
    it++;

    cudaEventElapsedTime(&ttmp, t[it], t[it+1]);
    s->time.transfer_to_host += ttmp;
    s->time.total2 += ttmp;
    it++;

    cudaEventElapsedTime(&ttmp, t[0], t[it+1]);
    s->time.total += ttmp;

    /* Cleanup */
    for (it=0; it<NT; it++) cudaEventDestroy(t[it]);

}


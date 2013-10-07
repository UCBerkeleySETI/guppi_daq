//# Copyright (C) 2013 Associated Universities, Inc. Washington DC, USA.
//# 
//# This program is free software; you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or
//# (at your option) any later version.
//# 
//# This program is distributed in the hope that it will be useful, but
//# WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//# General Public License for more details.
//# 
//# You should have received a copy of the GNU General Public License
//# along with this program; if not, write to the Free Software
//# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//# 
//# Correspondence concerning GBT software should be addressed as follows:
//#	GBT Operations
//#	National Radio Astronomy Observatory
//#	P. O. Box 2
//#	Green Bank, WV 24944-0002 USA

static char rcs_id[] =  "$Id$";

// Parent
#include "sse_transpose.h"
// Other
#ifdef __SSE__
#include "xmmintrin.h"
#endif//__SSE__

#include <stdio.h>

#define restrict __restrict__

#ifdef __SSE__ 
#define REGS __m128 row0,row1,row2,row3 
#else
#define REGS 
#endif

#define SSE_TILE_SIZE (4)

#define UNALIGNED_LOAD \
                row0 = _mm_loadu_ps(&in[inidx+xdim*0]);\
                row1 = _mm_loadu_ps(&in[inidx+xdim*1]);\
                row2 = _mm_loadu_ps(&in[inidx+xdim*2]);\
                row3 = _mm_loadu_ps(&in[inidx+xdim*3])
    
#define ALIGNED_LOAD \
                row0 = _mm_load_ps(&in[inidx+xdim*0]);\
                row1 = _mm_load_ps(&in[inidx+xdim*1]);\
                row2 = _mm_load_ps(&in[inidx+xdim*2]);\
                row3 = _mm_load_ps(&in[inidx+xdim*3])

#define UNALIGNED_STORE \
                _mm_storeu_ps(&out[outidx+ydim*0], row0);\
                _mm_storeu_ps(&out[outidx+ydim*1], row1);\
                _mm_storeu_ps(&out[outidx+ydim*2], row2);\
                _mm_storeu_ps(&out[outidx+ydim*3], row3)

#define ALIGNED_STORE \
                _mm_store_ps(&out[outidx+ydim*0], row0);\
                _mm_store_ps(&out[outidx+ydim*1], row1);\
                _mm_store_ps(&out[outidx+ydim*2], row2);\
                _mm_store_ps(&out[outidx+ydim*3], row3)

#define BEGIN_LOOP \
    do \
    { \
        int tile_x, tile_y; \
    for(tile_x = 0; tile_x < n_xtiles; ++tile_x)\
    {\
        for(tile_y = 0; tile_y < n_ytiles; ++tile_y)\
        {\
            REGS; \
                int dx = tile_x*SSE_TILE_SIZE; \
                int dy = tile_y*SSE_TILE_SIZE; \
                int inidx =  dx + dy*xdim; \
                int outidx = dy + dx*ydim;

#define END_LOOP \
        }\
    } \
    } while(0)

/**
 * Transpose a 4x4 block of float values. 
 *
 * The x86 architecture has super-scalar extensions (SSE) available
 * for processing vector data. The instructions have optimizations
 * which require 16 byte data alignment. 
 *
 * This function handles all of the alignment cases, using the most
 * efficient instructions for the alignment. If SSE instructions are
 * not available, non-SSE instructions will be used. 
 * While the alignment cases are handled at runtime, no effort
 * is made here to detect availablilty of SSE at run-time. (Just about any recent CPU
 * has the required SSE instructions available.)
 *
 * Note: out_ydim will be 0 for normal transposes. A non-zero value redefines ydim
 *       for application specific use. It also must be either zero or an integral
 *       multiple of SSE_TILE_SIZE.
 * Note: The routine will not work for dimensions which are not an integral multiple
 *       of the SSE_TILE_SIZE (4 for SSE instructions, 8 for AVX, not implemented here.) 
 *
 */
 


void sse_transpose(float * restrict in, float * restrict out, int xdim, int ydim, int out_ydim)
{

    enum Alignment { BOTH_ALIGNED, INPUT_UNALIGNED, OUTPUT_UNALIGNED, BOTH_UNALIGNED,
                     TILE_SIZE_ERROR };
                     
    int alignment_case=BOTH_ALIGNED;                     
    int n_xtiles = xdim/SSE_TILE_SIZE;
    int n_ytiles = ydim/SSE_TILE_SIZE;


    if ((long)in & 0xFL)
        alignment_case|=INPUT_UNALIGNED;
    if ((long)out & 0xFL)
        alignment_case|=OUTPUT_UNALIGNED;

    if (n_xtiles*SSE_TILE_SIZE != xdim || n_ytiles*SSE_TILE_SIZE != ydim ||
        (out_ydim != 0 && out_ydim%SSE_TILE_SIZE != 0))
        alignment_case = TILE_SIZE_ERROR; 
               
    /* This is an odd case where the output matrix is actually much larger
     * than the input. (i.e. the transpose is taking all of the input
     * and calculating a portion of the output. (application specific)
     */
     if (out_ydim != 0 && out_ydim%SSE_TILE_SIZE == 0)
         ydim = out_ydim;
         

#ifdef __SSE__
    switch (alignment_case)
    { 
    case BOTH_ALIGNED:
        BEGIN_LOOP;
        ALIGNED_LOAD;
        _MM_TRANSPOSE4_PS(row0, row1, row2, row3);
        ALIGNED_STORE;
        END_LOOP;
    break;
    case OUTPUT_UNALIGNED:
        BEGIN_LOOP;
        ALIGNED_LOAD;
        _MM_TRANSPOSE4_PS(row0, row1, row2, row3);
        UNALIGNED_STORE;
        END_LOOP;
    break;
    case INPUT_UNALIGNED:
        BEGIN_LOOP;
        UNALIGNED_LOAD;
        _MM_TRANSPOSE4_PS(row0, row1, row2, row3);
        ALIGNED_STORE;
        END_LOOP;
    break;
    case BOTH_UNALIGNED:
        BEGIN_LOOP;
        UNALIGNED_LOAD;
        _MM_TRANSPOSE4_PS(row0, row1, row2, row3);
        UNALIGNED_STORE;
        END_LOOP;
    break;
    case TILE_SIZE_ERROR:
        fprintf(stderr, "sse_transpose: dims must be a multiple of "
                        "block size (%d) xdim=%d ydim=%d\n",SSE_TILE_SIZE, xdim, ydim);
        return;
    break;
    }

#else
    BEGIN_LOOP;
            // If the SSE intrinsics aren't available,
            // do it the old fashioned way

            // -- XX XX XX //      // -- -- -- -- //
            // -- XX XX XX // ---\ // XX XX XX XX //
            // -- XX XX XX // ---/ // XX XX XX XX //
            // -- XX XX XX //      // XX XX XX XX //
            out[outidx+ydim*0+0] = in[inidx+xdim*0+0];
            out[outidx+ydim*0+1] = in[inidx+xdim*1+0];
            out[outidx+ydim*0+2] = in[inidx+xdim*2+0];
            out[outidx+ydim*0+3] = in[inidx+xdim*3+0];

            // XX -- XX XX //      // XX XX XX XX //
            // XX -- XX XX // ---\ // -- -- -- -- //
            // XX -- XX XX // ---/ // XX XX XX XX //
            // XX -- XX XX //      // XX XX XX XX //
            out[outidx+ydim*1+0] = in[inidx+xdim*0+1];
            out[outidx+ydim*1+1] = in[inidx+xdim*1+1];
            out[outidx+ydim*1+2] = in[inidx+xdim*2+1];
            out[outidx+ydim*1+3] = in[inidx+xdim*3+1];

            // XX XX -- XX //      // XX XX XX XX //
            // XX XX -- XX // ---\ // XX XX XX XX //
            // XX XX -- XX // ---/ // -- -- -- -- //
            // XX XX -- XX //      // XX XX XX XX //
            out[outidx+ydim*2+0] = in[inidx+xdim*0+2];
            out[outidx+ydim*2+1] = in[inidx+xdim*1+2];
            out[outidx+ydim*2+2] = in[inidx+xdim*2+2];
            out[outidx+ydim*2+3] = in[inidx+xdim*3+2];

            // XX XX XX -- //      // XX XX XX XX //
            // XX XX XX -- // ---\ // XX XX XX XX //
            // XX XX XX -- // ---/ // XX XX XX XX //
            // XX XX XX -- //      // -- -- -- -- //
            out[outidx+ydim*3+0] = in[inidx+xdim*0+3];
            out[outidx+ydim*3+1] = in[inidx+xdim*1+3];
            out[outidx+ydim*3+2] = in[inidx+xdim*2+3];
            out[outidx+ydim*3+3] = in[inidx+xdim*3+3];
    END_LOOP;
#endif//__SSE__
}

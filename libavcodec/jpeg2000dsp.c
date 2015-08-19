/*
 * JPEG 2000 DSP functions
 * Copyright (c) 2007 Kamil Nowosad
 * Copyright (c) 2013 Nicolas Bertrand <nicoinattendu@gmail.com>
 *
 * This file is part of Libav.
 *
 * Libav is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * Libav is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Libav; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "config.h"
#include "libavutil/attributes.h"
#include "jpeg2000dsp.h"

#ifdef __SSE__
#include <xmmintrin.h>
#endif

/* Inverse ICT parameters in float and integer.
 * int value = (float value) * (1<<16) */
static const float f_ict_params[4] = {
    1.402f,
    0.34413f,
    0.71414f,
    1.772f
};

static const int i_ict_params[4] = {
     91881,
     22553,
     46802,
    116130
};

static void mct_decode_sse(
                float* restrict c0,
                float* restrict c1,
                float* restrict c2,
                int n)
{
    int i;
    __m128 vrv, vgu, vgv, vbu;
    vrv = _mm_set1_ps(1.402f);
    vgu = _mm_set1_ps(0.34413f);
    vgv = _mm_set1_ps(0.71414f);
    vbu = _mm_set1_ps(1.772f);
    for (i = 0; i < (n >> 3); ++i) {
        __m128 vy, vu, vv;
        __m128 vr, vg, vb;

        vy = _mm_load_ps(c0);
        vu = _mm_load_ps(c1);
        vv = _mm_load_ps(c2);
        vr = _mm_add_ps(vy, _mm_mul_ps(vv, vrv));
        vg = _mm_sub_ps(_mm_sub_ps(vy, _mm_mul_ps(vu, vgu)), _mm_mul_ps(vv, vgv));
        vb = _mm_add_ps(vy, _mm_mul_ps(vu, vbu));
        _mm_store_ps(c0, vr);
        _mm_store_ps(c1, vg);
        _mm_store_ps(c2, vb);
        c0 += 4;
        c1 += 4;
        c2 += 4;
        vy = _mm_load_ps(c0);
        vu = _mm_load_ps(c1);
        vv = _mm_load_ps(c2);
        vr = _mm_add_ps(vy, _mm_mul_ps(vv, vrv));
        vg = _mm_sub_ps(_mm_sub_ps(vy, _mm_mul_ps(vu, vgu)), _mm_mul_ps(vv, vgv));
        vb = _mm_add_ps(vy, _mm_mul_ps(vu, vbu));
        _mm_store_ps(c0, vr);
        _mm_store_ps(c1, vg);
        _mm_store_ps(c2, vb);
        c0 += 4;
        c1 += 4;
        c2 += 4;
    }
    n &= 7;
}

static void ict_float(void *_src0, void *_src1, void *_src2, int csize)
{
   float *src0 = _src0, *src1 = _src1, *src2 = _src2;
#ifdef __SSE__
    mct_decode_sse(src0, src1, src2, csize);
#else
    float i0f, i1f, i2f;
    int i;

    for (i = 0; i < csize; i++) {
        i0f = *src0 + (f_ict_params[0] * *src2);
        i1f = *src0 - (f_ict_params[1] * *src1)
                    - (f_ict_params[2] * *src2);
        i2f = *src0 + (f_ict_params[3] * *src1);
        *src0++ = i0f;
        *src1++ = i1f;
        *src2++ = i2f;
    }
#endif
}

static void ict_int(void *_src0, void *_src1, void *_src2, int csize)
{
    int32_t *src0 = _src0, *src1 = _src1, *src2 = _src2;
    int32_t i0, i1, i2;
    int i;

    for (i = 0; i < csize; i++) {
        i0 = *src0 + (((i_ict_params[0] * *src2) + (1 << 15)) >> 16);
        i1 = *src0 - (((i_ict_params[1] * *src1) + (1 << 15)) >> 16)
                   - (((i_ict_params[2] * *src2) + (1 << 15)) >> 16);
        i2 = *src0 + (((i_ict_params[3] * *src1) + (1 << 15)) >> 16);
        *src0++ = i0;
        *src1++ = i1;
        *src2++ = i2;
    }
}

static void rct_int(void *_src0, void *_src1, void *_src2, int csize)
{
    int32_t *src0 = _src0, *src1 = _src1, *src2 = _src2;
    int32_t i0, i1, i2;
    int i;

    for (i = 0; i < csize; i++) {
        i1 = *src0 - (*src2 + *src1 >> 2);
        i0 = i1 + *src2;
        i2 = i1 + *src1;
        *src0++ = i0;
        *src1++ = i1;
        *src2++ = i2;
    }
}

av_cold void ff_jpeg2000dsp_init(Jpeg2000DSPContext *c)
{
    c->mct_decode[FF_DWT97]     = ict_float;
    c->mct_decode[FF_DWT53]     = rct_int;
    c->mct_decode[FF_DWT97_INT] = ict_int;
}

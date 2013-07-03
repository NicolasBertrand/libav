/*
 * Discrete wavelet transform
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

/**
 * @file
 * Discrete wavelet transform
 */
#ifdef __SSE__
#include <xmmintrin.h>
#endif

#include "libavutil/common.h"
#include "libavutil/mem.h"
#include "jpeg2000dwt.h"
#include "internal.h"

/* Defines for 9/7 DWT lifting parameters.
 * Parameters are in float. */
#define F_LFTG_ALPHA  1.586134342059924f
#define F_LFTG_BETA   0.052980118572961f
#define F_LFTG_GAMMA  0.882911075530934f
#define F_LFTG_DELTA  0.443506852043971f
#define F_LFTG_K      1.230174104914001f
#define F_LFTG_X      1.625732422f
/* FIXME: Why use 1.625732422 instead of 1/F_LFTG_K?
 * Incorrect value in JPEG2000 norm.
 * see (ISO/IEC 15444:1 (version 2002) F.3.8.2 */

/* Lifting parameters in integer format.
 * Computed as param = (float param) * (1 << 16) */
#define I_LFTG_ALPHA  103949
#define I_LFTG_BETA     3472
#define I_LFTG_GAMMA   57862
#define I_LFTG_DELTA   29066
#define I_LFTG_K       80621
#define I_LFTG_X      106544

static const float opj_dwt_alpha =  1.586134342f; /*  12994 */
static const float opj_dwt_beta  =  0.052980118f; /*    434 */
static const float opj_dwt_gamma = -0.882911075f; /*  -7233 */
static const float opj_dwt_delta = -0.443506852f; /*  -3633 */

static const float opj_K      = 1.230174105f; /*  10078 */
static const float opj_c13318 = 1.625732422f;



typedef union {
    float   f[4];
} opj_v4_t;

typedef struct v4dwt_local {
    opj_v4_t*   wavelet ;
    int32_t     dn ;
    int32_t     sn ;
    int32_t     cas ;
} opj_v4dwt_t ;
static void opj_v4dwt_decode(opj_v4dwt_t* restrict dwt);

static void opj_v4dwt_interleave_h(opj_v4dwt_t* restrict w, float* restrict a, int32_t x, int32_t size);

static void opj_v4dwt_interleave_v(opj_v4dwt_t* restrict v , float* restrict a , int32_t x, int32_t nb_elts_read);

#ifdef __SSE__
static void opj_v4dwt_decode_step1_sse(opj_v4_t* w, int32_t count, const __m128 c);

static void opj_v4dwt_decode_step2_sse(opj_v4_t* l, opj_v4_t* w, int32_t k, int32_t m, __m128 c);

#else
static void opj_v4dwt_decode_step1(opj_v4_t* w, int32_t count, const float c);

static void opj_v4dwt_decode_step2(opj_v4_t* l, opj_v4_t* w, int32_t k, int32_t m, float c);

#endif

static inline int32_t opj_int_min(int32_t a, int32_t b) {
    return a < b ? a : b;
}




static inline void extend53(int *p, int i0, int i1)
{
    p[i0 - 1] = p[i0 + 1];
    p[i1]     = p[i1 - 2];
    p[i0 - 2] = p[i0 + 2];
    p[i1 + 1] = p[i1 - 3];
}

static inline void extend97_float(float *p, int i0, int i1)
{
    int i;

    for (i = 1; i <= 4; i++) {
        p[i0 - i]     = p[i0 + i];
        p[i1 + i - 1] = p[i1 - i - 1];
    }
}

static inline void extend97_int(int32_t *p, int i0, int i1)
{
    int i;

    for (i = 1; i <= 4; i++) {
        p[i0 - i]     = p[i0 + i];
        p[i1 + i - 1] = p[i1 - i - 1];
    }
}

static void sr_1d53(int *p, int i0, int i1)
{
    int i;

    if (i1 == i0 + 1)
        return;

    extend53(p, i0, i1);

    for (i = i0 / 2; i < i1 / 2 + 1; i++)
        p[2 * i] -= (p[2 * i - 1] + p[2 * i + 1] + 2) >> 2;
    for (i = i0 / 2; i < i1 / 2; i++)
        p[2 * i + 1] += (p[2 * i] + p[2 * i + 2]) >> 1;
}

static void dwt_decode53(DWTContext *s, int *t)
{
    int lev;
    int w     = s->linelen[s->ndeclevels - 1][0];
    int32_t *line = s->i_linebuf;
    line += 3;

    for (lev = 0; lev < s->ndeclevels; lev++) {
        int lh = s->linelen[lev][0],
            lv = s->linelen[lev][1],
            mh = s->mod[lev][0],
            mv = s->mod[lev][1],
            lp;
        int *l;

        // HOR_SD
        l = line + mh;
        for (lp = 0; lp < lv; lp++) {
            int i, j = 0;
            // copy with interleaving
            for (i = mh; i < lh; i += 2, j++)
                l[i] = t[w * lp + j];
            for (i = 1 - mh; i < lh; i += 2, j++)
                l[i] = t[w * lp + j];

            sr_1d53(line, mh, mh + lh);

            for (i = 0; i < lh; i++)
                t[w * lp + i] = l[i];
        }

        // VER_SD
        l = line + mv;
        for (lp = 0; lp < lh; lp++) {
            int i, j = 0;
            // copy with interleaving
            for (i = mv; i < lv; i += 2, j++)
                l[i] = t[w * j + lp];
            for (i = 1 - mv; i < lv; i += 2, j++)
                l[i] = t[w * j + lp];

            sr_1d53(line, mv, mv + lv);

            for (i = 0; i < lv; i++)
                t[w * i + lp] = l[i];
        }
    }
}

static void sr_1d97_float(float *p, int i0, int i1)
{
    int i;

    if (i1 == i0 + 1)
        return;

    extend97_float(p, i0, i1);

    for (i = i0 / 2 - 1; i < i1 / 2 + 2; i++)
        p[2 * i]     -= F_LFTG_DELTA * (p[2 * i - 1] + p[2 * i + 1]);
    /* step 4 */
    for (i = i0 / 2 - 1; i < i1 / 2 + 1; i++)
        p[2 * i + 1] -= F_LFTG_GAMMA * (p[2 * i]     + p[2 * i + 2]);
    /*step 5*/
    for (i = i0 / 2; i < i1 / 2 + 1; i++)
        p[2 * i]     += F_LFTG_BETA  * (p[2 * i - 1] + p[2 * i + 1]);
    /* step 6 */
    for (i = i0 / 2; i < i1 / 2; i++)
        p[2 * i + 1] += F_LFTG_ALPHA * (p[2 * i]     + p[2 * i + 2]);
}

static void dwt_decode97_float(DWTContext *s, float *t)
{
    int lev;
    int w       = s->linelen[s->ndeclevels - 1][0];
    float *line = s->f_linebuf;
    float *data = t;
    /* position at index O of line range [0-5,w+5] cf. extend function */
    line += 5;

    for (lev = 0; lev < s->ndeclevels; lev++) {
        int lh = s->linelen[lev][0],
            lv = s->linelen[lev][1],
            mh = s->mod[lev][0],
            mv = s->mod[lev][1],
            lp;
        float *l;
        // HOR_SD
        l = line + mh;
        for (lp = 0; lp < lv; lp++) {
            int i, j = 0;
            // copy with interleaving
            for (i = mh; i < lh; i += 2, j++)
                l[i] = data[w * lp + j] * F_LFTG_K;
            for (i = 1 - mh; i < lh; i += 2, j++)
                l[i] = data[w * lp + j] * F_LFTG_X;

            sr_1d97_float(line, mh, mh + lh);

            for (i = 0; i < lh; i++)
                data[w * lp + i] = l[i];
        }

        // VER_SD
        l = line + mv;
        for (lp = 0; lp < lh; lp++) {
            int i, j = 0;
            // copy with interleaving
            for (i = mv; i < lv; i += 2, j++)
                l[i] = data[w * j + lp] * F_LFTG_K;
            for (i = 1 - mv; i < lv; i += 2, j++)
                l[i] = data[w * j + lp] * F_LFTG_X;

            sr_1d97_float(line, mv, mv + lv);

            for (i = 0; i < lv; i++)
                data[w * i + lp] = l[i];
        }
    }
}

static void opj_v4dwt_interleave_h(opj_v4dwt_t* restrict w, float* restrict a, int32_t x, int32_t size){
    float* restrict bi = (float*) (w->wavelet + w->cas);
    int count = w->sn;
    int i, k, idx;

    av_dlog(NULL, "interleave_h : count=%i bufsize=%i x=%i \n", count, size, x);
    for(k = 0; k < 2; ++k){
        if ( count + 3 * x < size && ((size_t) a & 0x0f) == 0 && ((size_t) bi & 0x0f) == 0 && (x & 0x0f) == 0 ) {
            /* Fast code path */
                        //av_dlog(NULL,"interleave_h: Fast code path\n");
            for(i = 0; i < count; ++i){
                int j = i;
                bi[i*8    ] = a[j];
                j += x;
                bi[i*8 + 1] = a[j];
                j += x;
                bi[i*8 + 2] = a[j];
                j += x;
                bi[i*8 + 3] = a[j];
            }
        }
        else {
            /* Slow code path */
                        //av_dlog(NULL,"interleave_h: Slow code path\n");
            for(i = 0; i < count; ++i){
                int j = i;
                bi[i*8    ] = a[j];
                j += x;
                if(j >= size) continue;
                bi[i*8 + 1] = a[j];
                j += x;
                if(j >= size) continue;
                bi[i*8 + 2] = a[j];
                j += x;
                if(j >= size) continue;
                bi[i*8 + 3] = a[j]; /* This one*/
            }
        }
        bi = (float*) (w->wavelet + 1 - w->cas);
        a += w->sn;
        size -= w->sn;
        count = w->dn;
    }
}

void opj_v4dwt_interleave_v(opj_v4dwt_t* restrict v , float* restrict a , int32_t x, int32_t nb_elts_read){
    opj_v4_t* restrict bi = v->wavelet + v->cas;
    int32_t i;

    for(i = 0; i < v->sn; ++i){
        memcpy(&bi[i*2], &a[i*x], nb_elts_read * sizeof(float));
    }

    a += v->sn * x;
    bi = v->wavelet + 1 - v->cas;

    for(i = 0; i < v->dn; ++i){
        memcpy(&bi[i*2], &a[i*x], nb_elts_read * sizeof(float));
    }
}

#ifdef __SSE__

void opj_v4dwt_decode_step1_sse(opj_v4_t* w, int32_t count, const __m128 c){
    __m128* restrict vw = (__m128*) w;
    int32_t i;
        //av_dlog(NULL,"STEP1 count = %i\n", count);
    /* 4x unrolled loop */
    for(i = 0; i < count >> 2; ++i){
        *vw = _mm_mul_ps(*vw, c);
        vw += 2;
        *vw = _mm_mul_ps(*vw, c);
        vw += 2;
        *vw = _mm_mul_ps(*vw, c);
        vw += 2;
        *vw = _mm_mul_ps(*vw, c);
        vw += 2;
    }
    count &= 3;
    for(i = 0; i < count; ++i){
        *vw = _mm_mul_ps(*vw, c);
        vw += 2;
    }
}

void opj_v4dwt_decode_step2_sse(opj_v4_t* l, opj_v4_t* w, int32_t k, int32_t m, __m128 c){
    __m128* restrict vl = (__m128*) l;
    __m128* restrict vw = (__m128*) w;
    int32_t i;
    __m128 tmp1, tmp2, tmp3;
    tmp1 = vl[0];
    for(i = 0; i < m; ++i){
        tmp2 = vw[-1];
        tmp3 = vw[ 0];
        vw[-1] = _mm_add_ps(tmp2, _mm_mul_ps(_mm_add_ps(tmp1, tmp3), c));
        tmp1 = tmp3;
        vw += 2;
    }
    vl = vw - 2;
    if(m >= k){
        return;
    }
    c = _mm_add_ps(c, c);
    c = _mm_mul_ps(c, vl[0]);
    for(; m < k; ++m){
        __m128 tmp = vw[-1];
        vw[-1] = _mm_add_ps(tmp, c);
        vw += 2;
    }
}

#else

void opj_v4dwt_decode_step1(opj_v4_t* w, int32_t count, const float c)
{
    float* restrict fw = (float*) w;
    int32_t i;
    for(i = 0; i < count; ++i){
        float tmp1 = fw[i*8    ];
        float tmp2 = fw[i*8 + 1];
        float tmp3 = fw[i*8 + 2];
        float tmp4 = fw[i*8 + 3];
        fw[i*8    ] = tmp1 * c;
        fw[i*8 + 1] = tmp2 * c;
        fw[i*8 + 2] = tmp3 * c;
        fw[i*8 + 3] = tmp4 * c;
    }
}

void opj_v4dwt_decode_step2(opj_v4_t* l, opj_v4_t* w, int32_t k, int32_t m, float c)
{
    float* restrict fl = (float*) l;
    float* restrict fw = (float*) w;
    int32_t i;
    for(i = 0; i < m; ++i){
        float tmp1_1 = fl[0];
        float tmp1_2 = fl[1];
        float tmp1_3 = fl[2];
        float tmp1_4 = fl[3];
        float tmp2_1 = fw[-4];
        float tmp2_2 = fw[-3];
        float tmp2_3 = fw[-2];
        float tmp2_4 = fw[-1];
        float tmp3_1 = fw[0];
        float tmp3_2 = fw[1];
        float tmp3_3 = fw[2];
        float tmp3_4 = fw[3];
        fw[-4] = tmp2_1 + ((tmp1_1 + tmp3_1) * c);
        fw[-3] = tmp2_2 + ((tmp1_2 + tmp3_2) * c);
        fw[-2] = tmp2_3 + ((tmp1_3 + tmp3_3) * c);
        fw[-1] = tmp2_4 + ((tmp1_4 + tmp3_4) * c);
        fl = fw;
        fw += 8;
    }
    if(m < k){
        float c1;
        float c2;
        float c3;
        float c4;
        c += c;
        c1 = fl[0] * c;
        c2 = fl[1] * c;
        c3 = fl[2] * c;
        c4 = fl[3] * c;
        for(; m < k; ++m){
            float tmp1 = fw[-4];
            float tmp2 = fw[-3];
            float tmp3 = fw[-2];
            float tmp4 = fw[-1];
            fw[-4] = tmp1 + c1;
            fw[-3] = tmp2 + c2;
            fw[-2] = tmp3 + c3;
            fw[-1] = tmp4 + c4;
            fw += 8;
        }
    }
}

#endif

void opj_v4dwt_decode(opj_v4dwt_t* restrict dwt)
{
    int32_t a, b;
    if(dwt->cas == 0) {
        if(!((dwt->dn > 0) || (dwt->sn > 1))){
            return;
        }
        a = 0;
        b = 1;
    }else{
        if(!((dwt->sn > 0) || (dwt->dn > 1))) {
            return;
        }
        a = 1;
        b = 0;
    }
        //av_dlog(NULL,"dwt->sn=%i, dwt->dn=%i \n", dwt->sn, dwt->dn);


#ifdef __SSE__
        //av_dlog(NULL," Calc with SSE\n");
    opj_v4dwt_decode_step1_sse(dwt->wavelet+a, dwt->sn, _mm_set1_ps(opj_K));
    opj_v4dwt_decode_step1_sse(dwt->wavelet+b, dwt->dn, _mm_set1_ps(opj_c13318));
    opj_v4dwt_decode_step2_sse(dwt->wavelet+b, dwt->wavelet+a+1, dwt->sn, opj_int_min(dwt->sn, dwt->dn-a), _mm_set1_ps(opj_dwt_delta));
    opj_v4dwt_decode_step2_sse(dwt->wavelet+a, dwt->wavelet+b+1, dwt->dn, opj_int_min(dwt->dn, dwt->sn-b), _mm_set1_ps(opj_dwt_gamma));
    opj_v4dwt_decode_step2_sse(dwt->wavelet+b, dwt->wavelet+a+1, dwt->sn, opj_int_min(dwt->sn, dwt->dn-a), _mm_set1_ps(opj_dwt_beta));
    opj_v4dwt_decode_step2_sse(dwt->wavelet+a, dwt->wavelet+b+1, dwt->dn, opj_int_min(dwt->dn, dwt->sn-b), _mm_set1_ps(opj_dwt_alpha));
#else
        //av_dlog(NULL," Calc without SSE\n");
    opj_v4dwt_decode_step1(dwt->wavelet+a, dwt->sn, opj_K);
    opj_v4dwt_decode_step1(dwt->wavelet+b, dwt->dn, opj_c13318);
    opj_v4dwt_decode_step2(dwt->wavelet+b, dwt->wavelet+a+1, dwt->sn, opj_int_min(dwt->sn, dwt->dn-a), opj_dwt_delta);
    opj_v4dwt_decode_step2(dwt->wavelet+a, dwt->wavelet+b+1, dwt->dn, opj_int_min(dwt->dn, dwt->sn-b), opj_dwt_gamma);
    opj_v4dwt_decode_step2(dwt->wavelet+b, dwt->wavelet+a+1, dwt->sn, opj_int_min(dwt->sn, dwt->dn-a), opj_dwt_beta);
    opj_v4dwt_decode_step2(dwt->wavelet+a, dwt->wavelet+b+1, dwt->dn, opj_int_min(dwt->dn, dwt->sn-b), opj_dwt_alpha);
#endif
}


int opj_dwt_decode(DWTContext *s, Jpeg2000Component *tilec)
{
    int numres = s->ndeclevels +1;
    opj_v4dwt_t h;
    opj_v4dwt_t v;

    Jpeg2000ResLevel* res = tilec->reslevel;

    int32_t rw = res->coord[0][1] - res->coord[0][0];    /* width of the resolution level computed */
    int32_t rh = res->coord[1][1] - res->coord[1][0];    /* height of the resolution level computed */

    int32_t w = tilec->coord[0][1] - tilec->coord[0][0];
    //TODO: temp usage of w
    h.wavelet = (opj_v4_t*) av_malloc((w+5) * sizeof(opj_v4_t));
    v.wavelet = h.wavelet;

    while( --numres) {
        float * restrict aj =  tilec->f_data;
        int32_t bufsize = (tilec->coord[0][1] - tilec->coord[0][0]) * (tilec->coord[1][1] - tilec->coord[1][0]);
        int32_t j;

        h.sn = rw;
        v.sn = rh;

        ++res;

        rw = res->coord[0][1] - res->coord[0][0];    /* width of the resolution level computed */
        rh = res->coord[1][1] - res->coord[1][0];    /* height of the resolution level computed */
        h.dn = rw - h.sn;
        h.cas = res->coord[0][0] % 2;
        for(j = rh; j > 3; j -= 4) {
            int32_t k;
            opj_v4dwt_interleave_h(&h, aj, w, bufsize);
            opj_v4dwt_decode(&h);
            for(k = rw; --k >= 0;){
                aj[k    ] = h.wavelet[k].f[0];
                aj[k+w  ] = h.wavelet[k].f[1];
                aj[k+w*2] = h.wavelet[k].f[2];
                aj[k+w*3] = h.wavelet[k].f[3];
            }
            aj += w*4;
            bufsize -= w*4;
        }

        if (rh & 0x03) {
            int32_t k;
            j = rh & 0x03;
            opj_v4dwt_interleave_h(&h, aj, w, bufsize);
            opj_v4dwt_decode(&h);
            for(k = rw; --k >= 0;){
                switch(j) {
                    case 3: aj[k+w*2] = h.wavelet[k].f[2];
                    case 2: aj[k+w  ] = h.wavelet[k].f[1];
                    case 1: aj[k    ] = h.wavelet[k].f[0];
                }
            }
        }
        v.dn = rh - v.sn;
        v.cas = res->coord[1][0] % 2;

        aj = tilec->f_data;
        for(j = rw; j > 3; j -= 4){
            int32_t k;

            opj_v4dwt_interleave_v(&v, aj, w, 4);
            opj_v4dwt_decode(&v);

            for(k = 0; k < rh; ++k){
                memcpy(&aj[k*w], &v.wavelet[k], 4 * sizeof(float));
            }
            aj += 4;
        }

        if (rw & 0x03){
            int32_t k;

            j = rw & 0x03;

            opj_v4dwt_interleave_v(&v, aj, w, j);
            opj_v4dwt_decode(&v);

            for(k = 0; k < rh; ++k){
                memcpy(&aj[k*w], &v.wavelet[k], j * sizeof(float));
            }
        }

    }
    av_freep(&h.wavelet);
    return 0;
}



static void sr_1d97_int(int32_t *p, int i0, int i1)
{
    int i;

    if (i1 == i0 + 1)
        return;

    extend97_int(p, i0, i1);

    for (i = i0 / 2 - 1; i < i1 / 2 + 2; i++)
        p[2 * i]     -= (I_LFTG_DELTA * (p[2 * i - 1] + p[2 * i + 1]) + (1 << 15)) >> 16;
    /* step 4 */
    for (i = i0 / 2 - 1; i < i1 / 2 + 1; i++)
        p[2 * i + 1] -= (I_LFTG_GAMMA * (p[2 * i]     + p[2 * i + 2]) + (1 << 15)) >> 16;
    /*step 5*/
    for (i = i0 / 2; i < i1 / 2 + 1; i++)
        p[2 * i]     += (I_LFTG_BETA  * (p[2 * i - 1] + p[2 * i + 1]) + (1 << 15)) >> 16;
    /* step 6 */
    for (i = i0 / 2; i < i1 / 2; i++)
        p[2 * i + 1] += (I_LFTG_ALPHA * (p[2 * i]     + p[2 * i + 2]) + (1 << 15)) >> 16;
}

static void dwt_decode97_int(DWTContext *s, int32_t *t)
{
    int lev;
    int w       = s->linelen[s->ndeclevels - 1][0];
    int32_t *line = s->i_linebuf;
    int32_t *data = t;
    /* position at index O of line range [0-5,w+5] cf. extend function */
    line += 5;

    for (lev = 0; lev < s->ndeclevels; lev++) {
        int lh = s->linelen[lev][0],
            lv = s->linelen[lev][1],
            mh = s->mod[lev][0],
            mv = s->mod[lev][1],
            lp;
        int32_t *l;
        // HOR_SD
        l = line + mh;
        for (lp = 0; lp < lv; lp++) {
            int i, j = 0;
            // rescale with interleaving
            for (i = mh; i < lh; i += 2, j++)
                l[i] = ((data[w * lp + j] * I_LFTG_K) + (1 << 15)) >> 16;
            for (i = 1 - mh; i < lh; i += 2, j++)
                l[i] = ((data[w * lp + j] * I_LFTG_X) + (1 << 15)) >> 16;

            sr_1d97_int(line, mh, mh + lh);

            for (i = 0; i < lh; i++)
                data[w * lp + i] = l[i];
        }

        // VER_SD
        l = line + mv;
        for (lp = 0; lp < lh; lp++) {
            int i, j = 0;
            // rescale with interleaving
            for (i = mv; i < lv; i += 2, j++)
                l[i] = ((data[w * j + lp] * I_LFTG_K) + (1 << 15)) >> 16;
            for (i = 1 - mv; i < lv; i += 2, j++)
                l[i] = ((data[w * j + lp] * I_LFTG_X) + (1 << 15)) >> 16;

            sr_1d97_int(line, mv, mv + lv);

            for (i = 0; i < lv; i++)
                data[w * i + lp] = l[i];
        }
    }
}

int ff_jpeg2000_dwt_init(DWTContext *s, uint16_t border[2][2],
                         int decomp_levels, int type)
{
    int i, j, lev = decomp_levels, maxlen,
        b[2][2];

    s->ndeclevels = decomp_levels;
    s->type       = type;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            b[i][j] = border[i][j];

    maxlen = FFMAX(b[0][1] - b[0][0],
                   b[1][1] - b[1][0]);
    while (--lev >= 0)
        for (i = 0; i < 2; i++) {
            s->linelen[lev][i] = b[i][1] - b[i][0];
            s->mod[lev][i]     = b[i][0] & 1;
            for (j = 0; j < 2; j++)
                b[i][j] = (b[i][j] + 1) >> 1;
        }
    switch (type) {
    case FF_DWT97:
        s->f_linebuf = av_malloc((maxlen + 12) * sizeof(*s->f_linebuf));
        if (!s->f_linebuf)
            return AVERROR(ENOMEM);
        break;
     case FF_DWT97_INT:
        s->i_linebuf = av_malloc((maxlen + 12) * sizeof(*s->i_linebuf));
        if (!s->i_linebuf)
            return AVERROR(ENOMEM);
        break;
    case FF_DWT53:
        s->i_linebuf = av_malloc((maxlen +  6) * sizeof(*s->i_linebuf));
        if (!s->i_linebuf)
            return AVERROR(ENOMEM);
        break;
    default:
        return -1;
    }
    return 0;
}

int ff_dwt_decode(DWTContext *s, void *t)
{
    switch (s->type) {
    case FF_DWT97:
        dwt_decode97_float(s, t);
        break;
    case FF_DWT97_INT:
        dwt_decode97_int(s, t);
        break;
    case FF_DWT53:
        dwt_decode53(s, t);
        break;
    default:
        return -1;
    }
    return 0;
}

void ff_dwt_destroy(DWTContext *s)
{
    av_freep(&s->f_linebuf);
    av_freep(&s->i_linebuf);
}

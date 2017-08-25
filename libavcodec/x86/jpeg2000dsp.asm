;******************************************************************************
;* SIMD-optimized JPEG2000 DSP functions
;* Copyright (c) 2014 Nicolas Bertrand
;* Copyright (c) 2015 James Almer
;* Copyright (c) 2017 Maxime Taisant
;*
;* This file is part of FFmpeg.
;*
;* FFmpeg is free software; you can redistribute it and/or
;* modify it under the terms of the GNU Lesser General Public
;* License as published by the Free Software Foundation; either
;* version 2.1 of the License, or (at your option) any later version.
;*
;* FFmpeg is distributed in the hope that it will be useful,
;* but WITHOUT ANY WARRANTY; without even the implied warranty of
;* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;* Lesser General Public License for more details.
;*
;* You should have received a copy of the GNU Lesser General Public
;* License along with FFmpeg; if not, write to the Free Software
;* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
;******************************************************************************

%include "libavutil/x86/x86util.asm"

SECTION_RODATA 32

pf_ict0: times 8 dd 1.402
pf_ict1: times 8 dd 0.34413
pf_ict2: times 8 dd 0.71414
pf_ict3: times 8 dd 1.772

F_LFTG_K: dd 1.230174104914001
F_LFTG_X: dd 0.812893066115961

F_LFTG_ALPHA: times 8 dd 1.586134342059924
F_LFTG_BETA: times 8 dd 0.052980118572961
F_LFTG_GAMMA: times 8 dd 0.882911075530934
F_LFTG_DELTA: times 8 dd 0.443506852043971

TWO: dd 2.0

SECTION .text

;***********************************************************************
; ff_ict_float_<opt>(float *src0, float *src1, float *src2, int csize)
;***********************************************************************
%macro ICT_FLOAT 1
cglobal ict_float, 4, 4, %1, src0, src1, src2, csize
    shl  csized, 2
    add   src0q, csizeq
    add   src1q, csizeq
    add   src2q, csizeq
    neg  csizeq
    movaps   m6, [pf_ict0]
    movaps   m7, [pf_ict1]
    %define ICT0 m6
    %define ICT1 m7

%if ARCH_X86_64
    movaps   m8, [pf_ict2]
    %define ICT2 m8
%if cpuflag(avx)
    movaps   m3, [pf_ict3]
    %define ICT3 m3
%else
    movaps   m9, [pf_ict3]
    %define ICT3 m9
%endif

%else ; ARCH_X86_32
    %define ICT2 [pf_ict2]
%if cpuflag(avx)
    movaps   m3, [pf_ict3]
    %define ICT3 m3
%else
    %define ICT3 [pf_ict3]
%endif

%endif ; ARCH

align 16
.loop:
    movaps   m0, [src0q+csizeq]
    movaps   m1, [src1q+csizeq]
    movaps   m2, [src2q+csizeq]

%if cpuflag(avx)
    mulps    m5, m1, ICT1
    mulps    m4, m2, ICT0
    mulps    m1, m1, ICT3
    mulps    m2, m2, ICT2
    subps    m5, m0, m5
%else ; sse
    movaps   m3, m1
    movaps   m4, m2
    movaps   m5, m0
    mulps    m3, ICT1
    mulps    m4, ICT0
    mulps    m1, ICT3
    mulps    m2, ICT2
    subps    m5, m3
%endif
    addps    m4, m4, m0
    addps    m0, m0, m1
    subps    m5, m5, m2

    movaps   [src0q+csizeq], m4
    movaps   [src2q+csizeq], m0
    movaps   [src1q+csizeq], m5
    add  csizeq, mmsize
    jl .loop
    REP_RET
%endmacro

INIT_XMM sse
ICT_FLOAT 10
INIT_YMM avx
ICT_FLOAT 9

;***************************************************************************
; ff_rct_int_<opt>(int32_t *src0, int32_t *src1, int32_t *src2, int csize)
;***************************************************************************
%macro RCT_INT 0
cglobal rct_int, 4, 4, 4, src0, src1, src2, csize
    shl  csized, 2
    add   src0q, csizeq
    add   src1q, csizeq
    add   src2q, csizeq
    neg  csizeq

align 16
.loop:
    mova   m1, [src1q+csizeq]
    mova   m2, [src2q+csizeq]
    mova   m0, [src0q+csizeq]
    paddd  m3, m1, m2
    psrad  m3, 2
    psubd  m0, m3
    paddd  m1, m0
    paddd  m2, m0
    mova   [src1q+csizeq], m0
    mova   [src2q+csizeq], m1
    mova   [src0q+csizeq], m2
    add  csizeq, mmsize
    jl .loop
    REP_RET
%endmacro

INIT_XMM sse2
RCT_INT
%if HAVE_AVX2_EXTERNAL
INIT_YMM avx2
RCT_INT
%endif

;***********************************************************************
; ff_sr_ld97_float_<opt>(float *line, int i0, int i1)
;***********************************************************************
%macro SR1D97FLOAT 1
cglobal sr_1d97_float, 3, 5, %1, line, i0, i1, j0, j1
    mov   j0q, i0q
    mov   j1q, i1q
    add   j0q, 1
    cmp   j1q, j0q
    jg %%extend
    sub   j0q, 2
    jnz %%else
    movss  m0, [lineq+4]
    movss  m1, [F_LFTG_K]
    movss  m2, [TWO]
    divss  m1, m2
    mulss  m0, m1
    movss  [lineq+4], m0
    jmp %%end

%%else:
    movss  m0, [lineq]
    movss  m1, [F_LFTG_X]
    mulss  m0, m1
    movss [lineq], m0
    jmp %%end

%%extend:
    shl   i0d, 2
    shl   i1d, 2
    mov   j0q, i0q
    mov   j1q, i1q
    movups m0, [lineq+j0q+4]
    shufps m0, m0, q0123
    movups [lineq+j0q-16], m0
    movups m0, [lineq+j1q-20]
    shufps m0, m0, q0123
    movups [lineq+j1q], m0

    movups m3, [F_LFTG_DELTA]
    mov   j0q, i0q
    mov   j1q, i1q
    shr   j0q, 1
    sub   j0q, 4
    shr   j1q, 1
    add   j1q, 8
    cmp   j0q, j1q
    jge %%beginloop2
%%loop1:
    add   j0q, 12
    cmp   j0q, j1q
    jge %%endloop1
 
    ;line{2*i,2*(i+1),2*(i+2),2*(i+3)} -= F_LFTG_DELTA*(line{2*i-1,2*(i+1)-1,2*(i+2)-1,2*(i+3)-1}+line{2*i+1,2*(i+1)+1,2*(i+2)+1,2*(i+3)+1})
    movups m0, [lineq+2*j0q-28]
    movups m4, [lineq+2*j0q-12]
    movaps m1, m0
    shufps m0, m4, q3131
    shufps m1, m4, q2020
    movups m2, [lineq+2*j0q-24]
    movups m5, [lineq+2*j0q-8] 
    shufps m2, m5, q3131
    addps  m2, m1
    mulps  m2, m3
    subps  m0, m2
    movaps m4, m1
    movlhps m1, m0
    shufps m1, m1, q3120
    shufps m4, m0, q3232
    shufps m4, m4, q3120
    movups [lineq+2*j0q-28], m1
    movups [lineq+2*j0q-12], m4

    add   j0q, 4
    cmp   j0q, j1q
    jge %%beginloop2
    jmp %%loop1
  
%%endloop1:
    sub   j0q, 12
%%littleloop1:
    movss  m0, [lineq+2*j0q]
    movss  m1, [lineq+2*j0q-4]
    movss  m2, [lineq+2*j0q+4]
    addss  m1, m2
    mulss  m1, m3
    subss  m0, m1
    movss [lineq+2*j0q], m0
    add   j0q, 4
    cmp   j0q, j1q
    jl %%littleloop1

%%beginloop2:
    movups m3, [F_LFTG_GAMMA]
    mov   j0q, i0q
    mov   j1q, i1q
    shr   j0q, 1
    sub   j0q, 4
    shr   j1q, 1
    add   j1q, 4
    cmp   j0q, j1q
    jge %%beginloop3
%%loop2:
    add   j0q, 12
    cmp   j0q, j1q
    jge %%endloop2
 
    ;line{2*i+1,2*(i+1)+1,2*(i+2)+1,2*(i+3)+1} -= F_LFTG_GAMMA*(line{2*i,2*(i+1),2*(i+2),2*(i+3)}+line{2*i+2,2*(i+1)+2,2*(i+2)+2,2*(i+3)+2})
    movups m0, [lineq+2*j0q-24]
    movups m4, [lineq+2*j0q-8]
    movaps m1, m0
    shufps m0, m4, q3131
    shufps m1, m4, q2020
    movups m2, [lineq+2*j0q-20]
    movups m5, [lineq+2*j0q-4] 
    shufps m2, m5, q3131
    addps  m2, m1
    mulps  m2, m3
    subps  m0, m2
    movaps m4, m1
    movlhps m1, m0
    shufps m1, m1, q3120
    shufps m4, m0, q3232
    shufps m4, m4, q3120
    movups [lineq+2*j0q-24], m1
    movups [lineq+2*j0q-8], m4

    add   j0q, 4
    cmp   j0q, j1q
    jge %%beginloop3
    jmp %%loop2
  
%%endloop2:
    sub   j0q, 12
%%littleloop2:
    movss  m0, [lineq+2*j0q+4]
    movss  m1, [lineq+2*j0q]
    movss  m2, [lineq+2*j0q+8]
    addss  m1, m2
    mulss  m1, m3
    subss  m0, m1
    movss  [lineq+2*j0q+4], m0
    add   j0q, 4
    cmp   j0q, j1q
    jl %%littleloop2

%%beginloop3:
    movups m3, [F_LFTG_BETA]
    mov   j0q, i0q
    mov   j1q, i1q
    shr   j0q, 1
    sub   j0q, 4
    shr   j1q, 1
    add   j1q, 8
    cmp   j0q, j1q
    jge %%beginloop4
%%loop3:
    add   j0q, 12
    cmp   j0q, j1q
    jge %%endloop3
 
    ;line{2*i,2*(i+1),2*(i+2),2*(i+3)} += F_LFTG_BETA*(line{2*i-1,2*(i+1)-1,2*(i+2)-1,2*(i+3)-1}+line{2*i+1,2*(i+1)+1,2*(i+2)+1,2*(i+3)+1})
    movups m0, [lineq+2*j0q-28]
    movups m4, [lineq+2*j0q-12]
    movaps m1, m0
    shufps m0, m4, q3131
    shufps m1, m4, q2020
    movups m2, [lineq+2*j0q-24]
    movups m5, [lineq+2*j0q-8] 
    shufps m2, m5, q3131
    addps  m2, m1
    mulps  m2, m3
    addps  m0, m2
    movaps m4, m1
    movlhps m1, m0
    shufps m1, m1, q3120
    shufps m4, m0, q3232
    shufps m4, m4, q3120
    movups [lineq+2*j0q-28], m1
    movups [lineq+2*j0q-12], m4

    add   j0q, 4
    cmp   j0q, j1q
    jge %%beginloop4
    jmp %%loop3
  
%%endloop3:
    sub   j0q, 12
%%littleloop3:
    movss  m0, [lineq+2*j0q]
    movss  m1, [lineq+2*j0q-4]
    movss  m2, [lineq+2*j0q+4]
    addss  m1, m2
    mulss  m1, m3
    addss  m0, m1
    movss [lineq+2*j0q], m0
    add   j0q, 4
    cmp   j0q, j1q
    jl %%littleloop3

%%beginloop4:
    movups m3, [F_LFTG_ALPHA]
    mov   j0q, i0q
    mov   j1q, i1q
    shr   j0q, 1
    sub   j0q, 4
    shr   j1q, 1
    add   j1q, 4
    cmp   j0q, j1q
    jge %%end
%%loop4:
    add   j0q, 12
    cmp   j0q, j1q
    jge %%endloop4
 
    ;line{2*i+1,2*(i+1)+1,2*(i+2)+1,2*(i+3)+1} += F_LFTG_ALPHA*(line{2*i,2*(i+1),2*(i+2),2*(i+3)}+line{2*i+2,2*(i+1)+2,2*(i+2)+2,2*(i+3)+2})
    movups m0, [lineq+2*j0q-24]
    movups m4, [lineq+2*j0q-8]
    movaps m1, m0
    shufps m0, m4, q3131
    shufps m1, m4, q2020
    movups m2, [lineq+2*j0q-20]
    movups m5, [lineq+2*j0q-4] 
    shufps m2, m5, q3131
    addps  m2, m1
    mulps  m2, m3
    addps  m0, m2
    movaps m4, m1
    movlhps m1, m0
    shufps m1, m1, q3120
    shufps m4, m0, q3232
    shufps m4, m4, q3120
    movups [lineq+2*j0q-24], m1
    movups [lineq+2*j0q-8], m4

    add   j0q, 4
    cmp   j0q, j1q
    jge %%end
    jmp %%loop4
  
%%endloop4:
    sub   j0q, 12
%%littleloop4:
    movss  m0, [lineq+2*j0q+4]
    movss  m1, [lineq+2*j0q]
    movss  m2, [lineq+2*j0q+8]
    addss  m1, m2
    mulss  m1, m3
    addss  m0, m1
    movss [lineq+2*j0q+4], m0
    add   j0q, 4
    cmp   j0q, j1q
    jl %%littleloop4

%%end:
    REP_RET
%endmacro
    
INIT_XMM sse
SR1D97FLOAT 6

%macro SR1D97FLOAT_ 5      ; p, i0, i1, tmp0, tmp1
    mov    %4, %2
    mov    %5, %3
    add    %4, 1
    cmp    %5, %4
    jg %%extend
    sub    %4, 2
    jnz %%else
    movss  m0, [%1+4]
    movss  m1, [F_LFTG_K]
    movss  m2, [TWO]
    divss  m1, m2
    mulss  m0, m1
    movss  [%1+4], m0
    jmp %%end

%%else:
    movss  m0, [%1]
    movss  m1, [F_LFTG_X]
    mulss  m0, m1
    movss  [%1], m0
    jmp %%end

%%extend:
    shl    %2, 2
    shl    %3, 2
    mov    %4, %2
    mov    %5, %3
    movups m0, [%1+%4+4]
    shufps m0, m0, q0123
    movups [%1+%4-16], m0
    movups m0, [%1+%5-20]
    shufps m0, m0, q0123
    movups [%1+%5], m0

    movups m3, [F_LFTG_DELTA]
    mov    %4, %2
    mov    %5, %3
    shr    %4, 1
    sub    %4, 4
    shr    %5, 1
    add    %5, 8
    cmp    %4, %5
    jge %%beginloop2
%%loop1:
    add    %4, 12
    cmp    %4, %5
    jge %%endloop1
 
    movups m0, [%1+2*%4-28]
    movups m4, [%1+2*%4-12]
    movaps m1, m0
    shufps m0, m4, q3131
    shufps m1, m4, q2020
    movups m2, [%1+2*%4-24]
    movups m5, [%1+2*%4-8] 
    shufps m2, m5, q3131
    addps  m2, m1
    mulps  m2, m3
    subps  m0, m2
    movaps m4, m1
    movlhps m1, m0
    shufps m1, m1, q3120
    shufps m4, m0, q3232
    shufps m4, m4, q3120
    movups [%1+2*%4-28], m1
    movups [%1+2*%4-12], m4

    add    %4, 4
    cmp    %4, %5
    jge %%beginloop2
    jmp %%loop1
  
%%endloop1:
    sub    %4, 12
%%littleloop1:
    movss  m0, [%1+2*%4]
    movss  m1, [%1+2*%4-4]
    movss  m2, [%1+2*%4+4]
    addss  m1, m2
    mulss  m1, m3
    subss  m0, m1
    movss  [%1+2*%4], m0
    add    %4, 4
    cmp    %4, %5
    jl %%littleloop1

%%beginloop2:
    movups m3, [F_LFTG_GAMMA]
    mov    %4, %2
    mov    %5, %3
    shr    %4, 1
    sub    %4, 4
    shr    %5, 1
    add    %5, 4
    cmp    %4, %5
    jge %%beginloop3
%%loop2:
    add    %4, 12
    cmp    %4, %5
    jge %%endloop2
 
    movups m0, [%1+2*%4-24]
    movups m4, [%1+2*%4-8]
    movaps m1, m0
    shufps m0, m4, q3131
    shufps m1, m4, q2020
    movups m2, [%1+2*%4-20]
    movups m5, [%1+2*%4-4] 
    shufps m2, m5, q3131
    addps  m2, m1
    mulps  m2, m3
    subps  m0, m2
    movaps m4, m1
    movlhps m1, m0
    shufps m1, m1, q3120
    shufps m4, m0, q3232
    shufps m4, m4, q3120
    movups [%1+2*%4-24], m1
    movups [%1+2*%4-8], m4

    add    %4, 4
    cmp    %4, %5
    jge %%beginloop3
    jmp %%loop2
  
%%endloop2:
    sub    %4, 12
%%littleloop2:
    movss  m0, [%1+2*%4+4]
    movss  m1, [%1+2*%4]
    movss  m2, [%1+2*%4+8]
    addss  m1, m2
    mulss  m1, m3
    subss  m0, m1
    movss  [%1+2*%4+4], m0
    add    %4, 4
    cmp    %4, %5
    jl %%littleloop2

%%beginloop3:
    movups m3, [F_LFTG_BETA]
    mov    %4, %2
    mov    %5, %3
    shr    %4, 1
    sub    %4, 4
    shr    %5, 1
    add    %5, 8
    cmp    %4, %5
    jge %%beginloop4
%%loop3:
    add    %4, 12
    cmp    %4, %5
    jge %%endloop3

    movups m0, [%1+2*%4-28]
    movups m4, [%1+2*%4-12]
    movaps m1, m0
    shufps m0, m4, q3131
    shufps m1, m4, q2020
    movups m2, [%1+2*%4-24]
    movups m5, [%1+2*%4-8] 
    shufps m2, m5, q3131
    addps  m2, m1
    mulps  m2, m3
    addps  m0, m2
    movaps m4, m1
    movlhps m1, m0
    shufps m1, m1, q3120
    shufps m4, m0, q3232
    shufps m4, m4, q3120
    movups [%1+2*%4-28], m1
    movups [%1+2*%4-12], m4

    add    %4, 4
    cmp    %4, %5
    jge %%beginloop4
    jmp %%loop3
  
%%endloop3:
    sub    %4, 12
%%littleloop3:
    movss  m0, [%1+2*%4]
    movss  m1, [%1+2*%4-4]
    movss  m2, [%1+2*%4+4]
    addss  m1, m2
    mulss  m1, m3
    addss  m0, m1
    movss  [%1+2*%4], m0
    add    %4, 4
    cmp    %4, %5
    jl %%littleloop3

%%beginloop4:
    movups m3, [F_LFTG_ALPHA]
    mov    %4, %2
    mov    %5, %3
    shr    %4, 1
    sub    %4, 4
    shr    %5, 1
    add    %5, 4
    cmp    %4, %5
    jge %%end
%%loop4:
    add    %4, 12
    cmp    %4, %5
    jge %%endloop4
 
    movups m0, [%1+2*%4-24]
    movups m4, [%1+2*%4-8]
    movaps m1, m0
    shufps m0, m4, q3131
    shufps m1, m4, q2020
    movups m2, [%1+2*%4-20]
    movups m5, [%1+2*%4-4] 
    shufps m2, m5, q3131
    addps  m2, m1
    mulps  m2, m3
    addps  m0, m2
    movaps m4, m1
    movlhps m1, m0
    shufps m1, m1, q3120
    shufps m4, m0, q3232
    shufps m4, m4, q3120
    movups [%1+2*%4-24], m1
    movups [%1+2*%4-8], m4

    add    %4, 4
    cmp    %4, %5
    jge %%end
    jmp %%loop4
  
%%endloop4:
    sub    %4, 12
%%littleloop4:
    movss  m0, [%1+2*%4+4]
    movss  m1, [%1+2*%4]
    movss  m2, [%1+2*%4+8]
    addss  m1, m2
    mulss  m1, m3
    addss  m0, m1
    movss [%1+2*%4+4], m0
    add    %4, 4
    cmp    %4, %5
    jl %%littleloop4

%%end:
    shr    %2, 2
    shr    %3, 2
%endmacro


;***********************************************************************
; ff_hor_sd_float_<opt>(float *line, float *data, int mh, int lh, int lv, int w)
;***********************************************************************
%macro HORSDFLOAT 1
cglobal hor_sd_float, 6, 12, %1, line, data, mh, lh, lv, w, l, lp, i0, i1, j0, j1
    mov    lq, mhq
    shl    lq, 2
    add    lq, lineq
    shl   lhq, 2
    
    mov   lpq, 0
%%mainloop:
    ;j0 = w*lp+j
    mov   j0q, wq
    imul  j0q, lpq

    ;j1 = (lh-mh+1)/2 + j0
    mov   j1q, lhq
    shr   j1q, 2
    sub   j1q, mhq
    add   j1q, 1
    shr   j1q, 1
    add   j1q, j0q

    shl   j0q, 2
    shl   j1q, 2

    ;i1 = 1-mh
    mov   i1q, 1
    sub   i1q, mhq
    shl   i1q, 2

    ;i0 = mh
    mov   i0q, mhq
    shl   i0q, 2
 
    cmp   i0q, i1q
    jg %%i1i0

;i0 < i1
    cmp   i1q, lhq
    jge %%i0
    
    add   i0q, 4
    cmp   i0q, i1q
    jne %%inci0
 
;i1 = i0+1
%%beginloopi0i1   
    sub   i0q, 4

%%loopi0i1:
    add   i1q, 24
    cmp   i1q, lhq
    jge %%endloopi0i1

    ;l{i0,i0+2,i0+4,i0+6} <- data[j0:j0+3]
    ;l{i0,i0+3,i0+5,i0+7} = l{i1,i1+2,i1+4,i1+6} <- data[j1:j1+3]
    movups m0, [dataq+j0q]
    movups m2, [dataq+j1q]
    movaps m1, m0
    movlhps m0, m2
    shufps m0, m0, q3120
    shufps m1, m2, q3232
    shufps m1, m1, q3120
    movups [lq+i0q], m0
    movups [lq+i0q+16], m1

    add   i1q, 8
    add   i0q, 32
    add   j0q, 16
    add   j1q, 16
    cmp   i1q, lhq
    jl %%loopi0i1  
    cmp   i0q, lhq
    jge %%sr_1d

    ;i1>=lh & i0<lh
    movss m0, [dataq+j0q]
    movss [lq+i0q], m0
    jmp %%sr_1d

;i1 + 6 >= lh
%%endloopi0i1:
    sub   i1q, 24
%%littleloopi0i1:

    ;l[i0] <- data[j0]
    ;l[i1] <- data[j1]
    movss m0, [dataq+j0q]
    movss m1, [dataq+j1q]
    movss [lq+i0q], m0
    movss [lq+i1q], m1

    add   i0q, 8
    add   i1q, 8
    add   j0q, 4
    add   j1q, 4
    cmp   i1q, lhq
    jl %%littleloopi0i1
    cmp   i0q, lhq
    jge %%sr_1d

    ;i1>=lh & i0<lh
    movss  m0, [dataq+j0q]
    movss  [lq+i0q], m0
    jmp %%sr_1d

;i1 < i0
%%i1i0:
    cmp   i0q, lhq
    jge %%i1
    
    add   i1q, 4
    cmp   i0q, i1q
    jne %%inci1

;i0 = i1+1
%%beginloopi1i0    
    sub   i1q, 4

%%loopi1i0:
    add   i0q, 24
    cmp   i0q, lhq
    jge %%endloopi1i0

    ;l{i1,i1+2,i1+4,i1+6} <- data[j1:j1+3]
    ;l{i1,i1+3,i1+5,i1+7} = l{i0,i0+2,i0+4,i0+6} <- data[j0:j0+3]
    movups m0, [dataq+j1q]
    movups m2, [dataq+j0q]
    movaps m1, m0
    movlhps m0, m2
    shufps m0, m0, q3120
    shufps m1, m2, q3232
    shufps m1, m1, q3120
    movups [lq+i1q], m0
    movups [lq+i1q+16], m1

    add   i0q, 8
    add   i1q, 32
    add   j0q, 16
    add   j1q, 16
    cmp   i0q, lhq
    jl %%loopi1i0  
    cmp   i1q, lhq
    jge %%sr_1d

    ;i0>=lh & i1<lh
    movss  m0, [dataq+j1q]
    movss  [lq+i1q], m0
    jmp %%sr_1d

%%endloopi1i0:
    sub   i1q, 24
%%littleloopi1i0:

    ;l[i0] <- data[j0]
    ;l[i1] <- data[j1]
    movss  m0, [dataq+j1q]
    movss  m1, [dataq+j0q]
    movss  [lq+i1q], m0
    movss  [lq+i0q], m1

    add   i0q, 8
    add   i1q, 8
    add   j0q, 4
    add   j1q, 4
    cmp   i1q, lhq
    jl %%littleloopi1i0
    cmp   i0q, lhq
    jge %%sr_1d

    ;i0>=lh & i1<lh
    movss  m0, [dataq+j0q]
    movss  [lq+i0q], m0
    jmp %%sr_1d

;i0<i1 & i1>=lh
%%i0:
    cmp   i0q, lhq
    jge %%sr_1d
    movss  m0, [dataq+j0q]
    movss  [lq+i0q], m0
    add   i0q, 8
    add   j0q, 4
    jmp %%i0

;i1<i0 & i0>=lh
%%i1:
    cmp   i1q, lhq
    jge %%sr_1d
    movss  m0, [dataq+j1q]
    movss [lq+i1q], m0
    add   i1q, 8
    add   j1q, 4
    jmp %%i1

;i0 < i1-1
%%inci0:
    cmp   i0q, lhq
    jge %%sr_1d
    movss  m0, [dataq+j0q]
    movss  [lq+i0q-4], m0
    add   i0q, 8
    add   j0q, 4
    cmp   i0q, i1q
    je %%beginloopi0i1
    jmp %%inci0

;i1 < i0-1
%%inci1:
    cmp   i1q, lhq
    jge %%sr_1d
    movss  m0, [dataq+j1q]
    movss  [lq+i1q-4], m0
    add   i1q, 8
    add   j1q, 4
    cmp   i0q, i1q
    je %%beginloopi1i0
    jmp %%inci1

%%sr_1d:
    mov   i0q, mhq
    mov   i1q, lhq
    shr   i1q, 2
    add   i1q, mhq
    SR1D97FLOAT_ lineq, i0q, i1q, j0q, j1q

    mov   i0q, 0
    cmp   i0q, lhq
    jge %%endmainloop
    mov   j0q, wq
    imul  j0q, lpq
    shl   j0q, 2
%%subloop3:
    add   i0q, 12
    cmp   i0q, lhq
    jge %%endsubloop3

    movups  m0, [lq+i0q-12]
    movups  [dataq+j0q], m0

    add   i0q, 4
    add   j0q, 16
    cmp   i0q, lhq
    jge %%endmainloop
    jmp %%subloop3  

%%endsubloop3:
    sub   i0q, 12
%%littlesubloop3:
    movss  m0, [lq+i0q]
    movss  [dataq+j0q], m0

    add   i0q, 4
    add   j0q, 4
    cmp   i0q, lhq
    jl %%littlesubloop3  

%%endmainloop:
    add   lpq, 1
    cmp   lpq, lvq
    jl %%mainloop

    REP_RET
%endmacro

INIT_XMM sse
HORSDFLOAT 6

;***********************************************************************
; ff_ver_sd_float_<opt>(float *line, float *data, int mv, int lh, int lv, int w)
;***********************************************************************
%macro VERSDFLOAT 1
cglobal ver_sd_float, 6, 12, %1, line, data, mv, lh, lv, w, lp, i0, i1, j0, j1, inc
    shl   mvq, 2
    add   lineq, mvq
    mov   incq, lvq
    add   incq, 12
    shl   incq, 2
    shl   lvq, 2
    shl   wq, 2    
    
    mov   lpq, 0

%%mainloop:
    ;j0 = w*j+lp
    mov   j0q, lpq

    add   lpq, 3
    cmp   lpq, lhq
    jge %%beginmainloop2

    shr   lvq, 2
    shr   wq, 2
    ;j1 = w*(lv-mv+1)/2 + j0
    mov   j1q, lvq
    sub   j1q, mvq
    add   j1q, 1
    shr   j1q, 1
    imul  j1q, wq
    add   j1q, j0q

    shl   lvq, 2
    shl   wq, 2
    shl   j1q, 2
    shl   j0q, 2

    ;i1 = 1-mv
    mov   i1q, 4
    sub   i1q, mvq

    ;i0 = mv
    mov   i0q, mvq
 
    cmp   i0q, i1q
    jg %%i1i0

;i0 < i1
    cmp   i1q, lvq
    jge %%i0
  
    add   i0q, 4
    cmp   i0q, i1q
    jne %%inci0
 
;i1 = i0+1
%%beginloopi0i1   
    sub   i0q, 4

%%loopi0i1:

    movss   m0, [dataq+j0q]
    movss   m1, [dataq+j1q]
    movss   [lineq+i0q], m0
    movss   [lineq+i0q+4], m1
    movss   m0, [dataq+j0q+4]
    movss   m1, [dataq+j1q+4]
    add     lineq, incq
    movss   [lineq+i0q], m0
    movss   [lineq+i0q+4], m1
    movss   m0, [dataq+j0q+8]
    movss   m1, [dataq+j1q+8]
    add     lineq, incq
    movss   [lineq+i0q], m0
    movss   [lineq+i0q+4], m1
    movss   m0, [dataq+j0q+12]
    movss   m1, [dataq+j1q+12]
    add     lineq, incq
    movss   [lineq+i0q], m0
    movss   [lineq+i0q+4], m1
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq

    add     i1q, 8
    add     i0q, 8
    add     j0q, wq
    add     j1q, wq
    cmp     i1q, lvq
    jl %%loopi0i1
    cmp     i0q, lvq
    jge %%sr_1d

    movss   m0, [dataq+j0q]
    movss   [lineq+i0q], m0
    add     lineq, incq
    movss   m0, [dataq+j0q+4]
    movss   [lineq+i0q], m0
    add     lineq, incq
    movss   m0, [dataq+j0q+8]
    movss   [lineq+i0q], m0
    add     lineq, incq
    movss   m0, [dataq+j0q+12]
    movss   [lineq+i0q], m0
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq
    jmp %%sr_1d

;i1 < i0
%%i1i0:
    cmp     i0q, lvq
    jge %%i1
    
    add     i1q, 4
    cmp     i0q, i1q
    jne %%inci1

;i0 = i1+1
%%beginloopi1i0    
    sub     i1q, 4

%%loopi1i0:

    movss   m0, [dataq+j1q]
    movss   m1, [dataq+j0q]
    movss   [lineq+i1q], m0
    movss   [lineq+i1q+4], m1
    movss   m0, [dataq+j1q+4]
    movss   m1, [dataq+j0q+4]
    add     lineq, incq
    movss   [lineq+i1q], m0
    movss   [lineq+i1q+4], m1
    movss   m0, [dataq+j1q+8]
    movss   m1, [dataq+j0q+8]
    add     lineq, incq
    movss   [lineq+i1q], m0
    movss   [lineq+i1q+4], m1
    movss   m0, [dataq+j1q+12]
    movss   m1, [dataq+j0q+12]
    add     lineq, incq
    movss   [lineq+i1q], m0
    movss   [lineq+i1q+4], m1
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq

    add     i1q, 8
    add     i0q, 8
    add     j0q, wq
    add     j1q, wq
    cmp     i0q, lvq
    jl %%loopi1i0
    cmp     i1q, lvq
    jge %%sr_1d

    movss   m0, [dataq+j1q]
    movss   [lineq+i1q], m0
    add     lineq, incq
    movss   m0, [dataq+j1q+4]
    movss   [lineq+i1q], m0
    add     lineq, incq
    movss   m0, [dataq+j1q+8]
    movss   [lineq+i1q], m0
    add     lineq, incq
    movss   m0, [dataq+j1q+12]
    movss   [lineq+i1q], m0
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq
    jmp %%sr_1d

;i0<i1 & i1>=lv
%%i0:
    cmp     i0q, lvq
    jge %%sr_1d
    movss   m0, [dataq+j0q]
    movss   [lineq+i0q], m0
    add     lineq, incq
    movss   m0, [dataq+j0q+4]
    movss   [lineq+i0q], m0
    add     lineq, incq
    movss   m0, [dataq+j0q+8]
    movss   [lineq+i0q], m0
    add     lineq, incq
    movss   m0, [dataq+j0q+12]
    movss   [lineq+i0q], m0
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq
    add     i0q, 8
    add     j0q, wq
    jmp %%i0

;i1<i0 & i0>=lh
%%i1:
    cmp     i1q, lvq
    jge %%sr_1d
    movss   m0, [dataq+j1q]
    movss   [lineq+i1q], m0
    add     lineq, incq
    movss   m0, [dataq+j1q+4]
    movss   [lineq+i1q], m0
    add     lineq, incq
    movss   m0, [dataq+j1q+8]
    movss   [lineq+i1q], m0
    add     lineq, incq
    movss   m0, [dataq+j1q+12]
    movss   [lineq+i1q], m0
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq
    add     i1q, 8
    add     j1q, wq
    jmp %%i1

;i0 < i1-1
%%inci0:
    cmp     i0q, lvq
    jge %%sr_1d
    movss   m0, [dataq+j0q]
    movss   [lineq+i0q-4], m0
    add     lineq, incq
    movss   m0, [dataq+j0q+4]
    movss   [lineq+i0q-4], m0
    add     lineq, incq
    movss   m0, [dataq+j0q+8]
    movss   [lineq+i0q-4], m0
    add     lineq, incq
    movss   m0, [dataq+j0q+12]
    movss   [lineq+i0q-4], m0
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq
    add     i0q, 8
    add     j0q, wq
    cmp     i0q, i1q
    je %%beginloopi0i1
    jmp %%inci0

;i1 < i0-1
%%inci1:
    cmp     i1q, lvq
    jge %%sr_1d
    movss   m0, [dataq+j1q]
    movss   [lineq+i1q-4], m0
    add     lineq, incq
    movss   m0, [dataq+j1q+4]
    movss   [lineq+i1q-4], m0
    add     lineq, incq
    movss   m0, [dataq+j1q+8]
    movss   [lineq+i1q-4], m0
    add     lineq, incq
    movss   m0, [dataq+j1q+12]
    movss   [lineq+i1q-4], m0
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq
    add     i1q, 8
    add     j1q, wq
    cmp     i0q, i1q
    je %%beginloopi1i0
    jmp %%inci1

%%sr_1d:
    sub     lineq, mvq
    mov     i0q, mvq
    mov     i1q, lvq
    add     i1q, mvq
    shr     i0q, 2
    shr     i1q, 2
    SR1D97FLOAT_ lineq, i0q, i1q, j0q, j1q
    add     lineq, incq
    mov     i0q, mvq
    mov     i1q, lvq
    add     i1q, mvq
    shr     i0q, 2
    shr     i1q, 2
    SR1D97FLOAT_ lineq, i0q, i1q, j0q, j1q
    add     lineq, incq
    mov     i0q, mvq
    mov     i1q, lvq
    add     i1q, mvq
    shr     i0q, 2
    shr     i1q, 2
    SR1D97FLOAT_ lineq, i0q, i1q, j0q, j1q
    add     lineq, incq
    mov     i0q, mvq
    mov     i1q, lvq
    add     i1q, mvq
    shr     i0q, 2
    shr     i1q, 2
    SR1D97FLOAT_ lineq, i0q, i1q, j0q, j1q
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq
    add     lineq, mvq

    mov     i0q, 0
    ;cmp     i0q, lvq
    ;jge %%endmainloop
    mov     j0q, lpq
    sub     j0q, 3
    shl     j0q, 2
%%loop3:
    add     i0q, 12
    cmp     i0q, lvq
    jge %%endloop3

    movups  m0, [lineq+i0q-12]
    add     lineq, incq
    movups  m1, [lineq+i0q-12]
    add     lineq, incq
    movups  m2, [lineq+i0q-12]
    add     lineq, incq
    movups  m3, [lineq+i0q-12]
    movaps  m4, m0
    movaps  m5, m2
    movlhps m0, m1
    movlhps m2, m3
    movaps  m6, m0
    shufps  m0, m2, q2020
    shufps  m6, m2, q3131
    movaps  m2, m6
    movhlps m1, m4
    movhlps m3, m5
    movaps  m6, m1
    shufps  m1, m3, q2020
    shufps  m6, m3, q3131
    movaps  m3, m6
    movups  [dataq+j0q], m0
    add     j0q, wq
    movups  [dataq+j0q], m2
    add     j0q, wq
    movups  [dataq+j0q], m1
    add     j0q, wq
    movups  [dataq+j0q], m3
    add     j0q, wq
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq

    add     i0q, 4
    cmp     i0q, lvq
    jge %%endmainloop
    jmp %%loop3  

%%endloop3:
    sub     i0q, 12

%%littleloop3:
    movss   m0, [lineq+i0q]
    movss   [dataq+j0q], m0
    add     lineq, incq
    movss   m0, [lineq+i0q]
    movss   [dataq+j0q+4], m0
    add     lineq, incq
    movss   m0, [lineq+i0q]
    movss   [dataq+j0q+8], m0
    add     lineq, incq
    movss   m0, [lineq+i0q]
    movss   [dataq+j0q+12], m0
    sub     lineq, incq
    sub     lineq, incq
    sub     lineq, incq

    add     i0q, 4
    add     j0q, wq
    cmp     i0q, lvq
    jl %%littleloop3  

%%endmainloop:
    add     lpq, 1
    cmp     lpq, lhq
    jl %%mainloop
    jmp %%end

%%beginmainloop2:
    sub     lpq, 3
%%mainloop2:
    ;j0 = w*j+lp
    mov   j0q, lpq
    shl   j0q, 2

    ;i0 = mv
    mov   i0q, mvq
 
    cmp     i0q, lvq
    jge %%beginloop5
%%loop4:
    movss   m0, [dataq+j0q]
    movss   [lineq+i0q], m0

    add     j0q, wq
    add     i0q, 8
    cmp     i0q, lvq
    jl      %%loop4

%%beginloop5:
    ;i0 = 1-mv
    mov   i0q, 4
    sub   i0q, mvq
    cmp     i0q, lvq
    jge %%sr_1d_2
%%loop5:
    movss   m0, [dataq+j0q]
    movss   [lineq+i0q], m0

    add     j0q, wq
    add     i0q, 8
    cmp     i0q, lvq
    jl      %%loop5

%%sr_1d_2:
    sub     lineq, mvq
    mov     i0q, mvq
    mov     i1q, lvq
    add     i1q, mvq
    shr     i1q, 2
    shr     i0q, 2
    SR1D97FLOAT_ lineq, i0q, i1q, j0q, j1q
    add     lineq, mvq

    mov     i0q, 0
    cmp     i0q, lvq
    jge %%endmainloop
    mov     j0q, lpq
    shl     j0q, 2
%%loop6:
    movss   m0, [lineq+i0q]
    movss   [dataq+j0q], m0

    add     j0q, wq
    add     i0q, 4
    cmp     i0q, lvq
    jl %%loop6 

%%endmainloop2:
    add   lpq, 1
    cmp   lpq, lhq
    jl %%mainloop2

%%end:
    REP_RET
%endmacro

INIT_XMM sse
VERSDFLOAT 6

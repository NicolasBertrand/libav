/*
 * Discrete wavelet transform
 * Copyright (c) 2007 Kamil Nowosad
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

#ifndef AVCODEC_JPEG2000DWT_H
#define AVCODEC_JPEG2000DWT_H

/**
 * @file
 * Discrete wavelet transform
 */

#include <stdint.h>
#include "jpeg2000.h"


/**
 * Initialize DWT.
 * @param s                 DWT context
 * @param border            coordinates of transformed region {{x0, x1}, {y0, y1}}
 * @param decomp_levels     number of decomposition levels
 * @param type              0 for DWT 9/7; 1 for DWT 5/3
 */
int ff_jpeg2000_dwt_init(DWTContext *s, uint16_t border[2][2],
                         int decomp_levels, int type);

int ff_dwt_decode(DWTContext *s, void *t);

void ff_dwt_destroy(DWTContext *s);

int opj_dwt_decode(DWTContext *s, Jpeg2000Component *tilec);
#endif /* AVCODEC_JPEG2000DWT_H */

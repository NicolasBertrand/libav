/*
 * XYZ to RGB filter
 * Copyright (c) 2012 Matthias Buercher
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
 * xyz2rgb filter
 * Converts from XYZ to RGB space
 * Useful to convert jpeg2000 files from MXF containers in DCP
 * The filter has no parameters
 */

#include <string.h>

#include "libavutil/pixdesc.h"
#include "libavutil/common.h"
#include "libavutil/internal.h"
#include "libavutil/intreadwrite.h"
#include "libavutil/imgutils.h"
#include "avfilter.h"
#include "formats.h"
#include "internal.h"
#include "video.h"

/* pre defined color-spaces gamma */
#define XYZ_GAMMA (2.6f)
#define RGB_GAMMA (2.2f)

typedef struct XYZ2RGBContext {
    int xyzgamma[4096];
    int rgbgamma[4096];
    int (* matrix)[3];
} XYZ2RGBContext;

/*
 * XYZ to RGB conversion matrix
 *  R            | 3.2404542  -1.5371385 -0.4985314 | X
 *  G   = 4095 * | -0.9692660  1.8760108  0.0415560 | Y
 *  B            | 0.0556434  -0.2040259  1.0572252 | Z
 *
 *  The matrix is stored in int values
 */
static int xyz2rgb_matrix[3][3] = {
    { 13270, -6295, -2041 },
    { -3969,  7682,   170 },
    {   228,  -835,  4329 }
};

static int query_formats(AVFilterContext *ctx)
{
    AVFilterFormats *formats;
    int ret = -1 ;

    if (ctx->inputs[0]) {
        formats = NULL;
        ret = ff_add_format(&formats, AV_PIX_FMT_XYZ12LE);
        ret = ff_add_format(&formats, AV_PIX_FMT_XYZ12BE);
        ff_formats_ref(formats, &ctx->inputs[0]->out_formats);
    }
    if (ctx->outputs[0]) {
        formats = NULL;
        ret = ff_add_format(&formats, AV_PIX_FMT_RGB48);
        ff_formats_ref(formats, &ctx->outputs[0]->in_formats);
    }
    return ret;
}

/*
 * The gamma values are precalculated in an array
 * XYZ uses projector gamma 2.6
 * sRGB uses gamma 2.2
 * The gamma function is the inverse power function, calculated in [0..1] and scaled to 12-bit fixed point
 *
 * The matrix multipliers are precalculated and scaled to 12bit bitdepth (4096 values)
 * For documentation see
 * http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
 * http://en.wikipedia.org/wiki/SRGB
*/
static int config_props(AVFilterLink *inlink)
{
    XYZ2RGBContext *settings = inlink->dst->priv;
    int i;
    double xyzgamma = XYZ_GAMMA;
    double rgbgamma = 1.0 / RGB_GAMMA;

    /* set gamma vectors */
    for (i = 0; i < 4096; i++) {
        settings->xyzgamma[i] =
            lrintf(pow(i / 4095.0, xyzgamma) * 4095.0);
        settings->rgbgamma[i] =
            lrintf(pow(i / 4095.0, rgbgamma) * 4095.0);
    }

    settings->matrix = xyz2rgb_matrix;
    return 0;
}

static int filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterLink   *outlink  = inlink->dst->outputs[0];
    XYZ2RGBContext *settings = inlink->dst->priv;

    AVFrame *out;
    uint8_t *inrow, *outrow;
    int i, j;
    int r, g, b, x, y, z;

    out = ff_get_video_buffer(outlink, outlink->w, outlink->h);
    if (!out) {
        av_frame_free(&in);
        return AVERROR(ENOMEM);
    }

    out->pts = in->pts;
    inrow  = in->data[0];
    outrow = out->data[0];

    /*
     * In case of XYZ (12-bit depth) format input
     * conversion to RGB24 (8-bit depth) --> loss of precision
     * FIXME: test/implement output RGB48 for display
     */
    for (i = 0; i < inlink->h; i++) {
        for (j = 0; j < inlink->w * 6; j += 6) {
            // read according endianness and scale from 16bit to 12bit
            if (inlink->format == AV_PIX_FMT_XYZ12LE) {
                x = AV_RL16(inrow + j)     >> 4;
                y = AV_RL16(inrow + j + 2) >> 4;
                z = AV_RL16(inrow + j + 4) >> 4;
            } else {
                x = AV_RB16(inrow + j)     >> 4;
                y = AV_RB16(inrow + j + 2) >> 4;
                z = AV_RB16(inrow + j + 4) >> 4;
           } 

            // convert from XYZ to XYZlinear
            x = settings->xyzgamma[x];
            y = settings->xyzgamma[y];
            z = settings->xyzgamma[z];

            // convert from XYZlinear to sRGBlinear
            r = settings->matrix[0][0] * x +
                settings->matrix[0][1] * y +
                settings->matrix[0][2] * z >> 12;
            g = settings->matrix[1][0] * x +
                settings->matrix[1][1] * y +
                settings->matrix[1][2] * z >> 12;
            b = settings->matrix[2][0] * x +
                settings->matrix[1][2] * y +
                settings->matrix[2][2] * z >> 12;

            // limit values to 12-bit depth
            r = av_clip_c(r,0,4095);
            g = av_clip_c(g,0,4095);
            b = av_clip_c(b,0,4095);

            // convert from sRGBlinear to RGB and scale from 12bit to 16bit
            r = settings->rgbgamma[r] << 4;
            g = settings->rgbgamma[g] << 4;
            b = settings->rgbgamma[b] << 4;

            // write according endianness
            if (inlink->format == AV_PIX_FMT_XYZ12LE) {
                AV_WL16(outrow + j,     r);
                AV_WL16(outrow + j + 2, g);
                AV_WL16(outrow + j + 4, b);
            } else {
                AV_WB16(outrow + j,     r);
                AV_WB16(outrow + j + 2, g);
                AV_WB16(outrow + j + 4, b);
            }
        }
    inrow += in->linesize[0];
    outrow += out->linesize[0];
    }
    av_frame_free(&in);
    return ff_filter_frame(outlink, out);
}

static const AVFilterPad avfilter_vf_xyz2rgb_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .filter_frame = filter_frame,
        .config_props = config_props,
        .min_perms    = AV_PERM_READ,
    },
    { NULL }
};

static const AVFilterPad avfilter_vf_xyz2rgb_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter avfilter_vf_xyz2rgb = {
    .name          = "xyz2rgb",
    .description   = NULL_IF_CONFIG_SMALL("Converts XYZ to RGB."),
    .priv_size     = sizeof(XYZ2RGBContext),
    .query_formats = query_formats,
    .inputs        = avfilter_vf_xyz2rgb_inputs,
    .outputs       = avfilter_vf_xyz2rgb_outputs
};

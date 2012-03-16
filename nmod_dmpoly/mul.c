/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <string.h>
#include "nmod_dmpoly.h"

#define FLAT_MAX_LENGTH 25
#define KS_MIN_SIZE 400
#define KS_MIN_DENSITY 0.1

static __inline__ void
__nmod_dmpoly_mul_recursive_mul(void * z, const void * x, long xlen,
                        const void * y, long ylen, int vars, nmod_t mod)
{
    if (vars == 1)
    {
        _nmod_poly_mul((mp_ptr) z, (mp_srcptr) x, xlen,
                                   (mp_srcptr) y, ylen, mod);
    }
    else
    {
        arr_srcptr xp, yp;
        arr_ptr zp;
        long i, j;

        xp = x;
        yp = y;
        zp = z;

        for (i = 0; i < xlen; i++)
            _nmod_dmpoly_mul(zp + i, xp + i, yp, vars - 1, mod);

        if (ylen != 1)
        {
            arr_struct tmp;
            _nmod_dmpoly_init(&tmp, vars - 1);

            for (i = 0; i < ylen - 1; i++)
                _nmod_dmpoly_mul(zp + xlen + i,
                    yp + i + 1, xp + xlen - 1, vars - 1, mod);

            for (i = 0; i < xlen - 1; i++)
            {
                for (j = 0; j < ylen - 1; j++)
                {
                    _nmod_dmpoly_mul(&tmp, yp + 1 + j,
                        xp + i, vars - 1, mod);
                    _nmod_dmpoly_add(zp + i + j + 1, zp + i + j + 1,
                        &tmp, vars - 1, mod);
                }
            }

            _nmod_dmpoly_clear(&tmp, vars - 1);
        }
    }
}

void
__nmod_dmpoly_mul(void * z, const void * x, long xlen,
                            const void * y, long ylen, int vars, nmod_t mod)
{
    long i, xterms, yterms, * xbounds, * ybounds, xmax, ymax;
    double xrect, yrect;

    if (vars == 1)
    {
        _nmod_poly_mul((mp_ptr) z, (mp_srcptr) x, xlen,
                                   (mp_srcptr) y, ylen, mod);
        return;
    }

    xbounds = flint_calloc(vars, sizeof(long));
    ybounds = flint_calloc(vars, sizeof(long));

    xterms = __nmod_dmpoly_rect_hull(xbounds, x, xlen, vars);
    yterms = __nmod_dmpoly_rect_hull(ybounds, y, ylen, vars);

    xrect = xmax = xbounds[0];
    for (i = 1; i < vars; i++)
    {
        xrect *= xbounds[i];
        xmax = FLINT_MAX(xmax, xbounds[i]);
    }

    yrect = ymax = ybounds[0];
    for (i = 1; i < vars; i++)
    {
        yrect *= ybounds[i];
        ymax = FLINT_MAX(ymax, ybounds[i]);
    }

    if ((xterms > xrect * KS_MIN_DENSITY) &&
        (yterms > yrect * KS_MIN_DENSITY) &&
        (xmax > KS_MIN_SIZE) &&
        (ymax > KS_MIN_SIZE))
    {
        __nmod_dmpoly_mul_KS(z, x, y, xbounds, ybounds, vars, mod);
    }
    else if (ymax < FLAT_MAX_LENGTH && xmax < FLAT_MAX_LENGTH && vars == 2)
    {
        __nmod_dmpoly_mul_flat(z, x, y, xbounds, ybounds, vars, mod);
    }
    else
    {
        __nmod_dmpoly_mul_recursive_mul(z, x, xlen, y, ylen, vars, mod);
    }

    flint_free(xbounds);
    flint_free(ybounds);
}

void
_nmod_dmpoly_mul(arr_ptr z, arr_srcptr x,
                                        arr_srcptr y, int vars, nmod_t mod)
{
    long xlen, ylen, zlen;

    xlen = x->length;
    ylen = y->length;

    if (xlen == 0 || ylen == 0)
    {
        _nmod_dmpoly_zero(z, vars);
        return;
    }

    if (z == x || z == y)
    {
        arr_struct tmp;
        _nmod_dmpoly_init(&tmp, vars);
        _nmod_dmpoly_mul(&tmp, x, y, vars, mod);
        _nmod_dmpoly_swap(z, &tmp);
        _nmod_dmpoly_clear(&tmp, vars);
        return;
    }

    zlen = xlen + ylen - 1;
    _nmod_dmpoly_fit_length(z, zlen, vars);
    if (xlen >= ylen)
        __nmod_dmpoly_mul(z->coeffs, x->coeffs, xlen,
                                            y->coeffs, ylen, vars, mod);
    else
        __nmod_dmpoly_mul(z->coeffs, y->coeffs, ylen,
                                            x->coeffs, xlen, vars, mod);
    z->length = zlen;
    _nmod_dmpoly_normalise(z, vars);
}

void
nmod_dmpoly_mul(nmod_dmpoly_t z, const nmod_dmpoly_t x,
                                    const nmod_dmpoly_t y)
{
    _nmod_dmpoly_mul(&z->arr, &x->arr, &y->arr, x->vars, x->mod);
}

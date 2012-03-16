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

#include "nmod_dmpoly.h"

void
__nmod_dmpoly_bit_pack(mp_ptr xp, const void * x, long xlen, int vars,
                            const long * stride, long bits)
{
    if (vars == 1)
    {
        _nmod_poly_bit_pack(xp, (mp_srcptr) x, xlen, bits);
    }
    else
    {
        long i;

        for (i = 0; i < xlen; i++)
        {
            __nmod_dmpoly_bit_pack(xp + stride[0] * i,
                ((arr_srcptr) x + i)->coeffs,
                ((arr_srcptr) x + i)->length, vars - 1, stride + 1, bits);
        }
    }
}

/* TODO: this currently allocates space for leading zeros, which is madness */
void
__nmod_dmpoly_bit_unpack(void * z, mp_srcptr zp, int vars,
    const long * zbounds, const long * stride, long bits, nmod_t mod)
{
    if (vars == 1)
    {
        _nmod_poly_bit_unpack((mp_ptr) z, zbounds[0], zp, bits, mod);
    }
    else
    {
        long i;
        arr_ptr x = (arr_ptr) z;

        for (i = 0; i < zbounds[0]; i++)
        {
            _nmod_dmpoly_fit_length(x + i, zbounds[1], vars - 1);
            __nmod_dmpoly_bit_unpack((x + i)->coeffs,
                zp + stride[0] * i, vars - 1, zbounds + 1, stride + 1, bits, mod);
            (x + i)->length = zbounds[1];
            _nmod_dmpoly_normalise(x + i, vars - 1);
        }
    }
}

void
__nmod_dmpoly_mul_KS(arr_ptr z, arr_srcptr x, arr_srcptr y,
    const long * xbounds, const long * ybounds, int vars, nmod_t mod)
{
    mp_ptr zp, xp, yp;
    long i, maxterms, bits;
    long * stride, * zbounds;
    long xlimbs, ylimbs;

    stride = flint_calloc(vars, sizeof(long));
    zbounds = flint_calloc(vars, sizeof(long));

    for (i = 0; i < vars; i++)
        zbounds[i] = (xbounds[i] + ybounds[i] - 1);

    /* bound coefficients.
       TODO: compute max bits instead of using modulus.
       TODO: can overflow if polynomial is extremely sparse (but then
             we shouldn't be using KS...) */
    maxterms = 1;
    for (i = 0; i < vars; i++)
        maxterms *= FLINT_MIN(xbounds[i], ybounds[i]);
    bits = FLINT_BIT_COUNT(maxterms) + 2*FLINT_BIT_COUNT(mod.n - 1);

    /* number of limbs at the bottom level */
    stride[vars-1] = (zbounds[vars-1] * bits + FLINT_BITS - 1) / FLINT_BITS;
    for (i = vars - 2; i >= 0; i--)
        stride[i] = stride[i + 1] * zbounds[i];

    xlimbs = stride[1] * xbounds[0];
    ylimbs = stride[1] * ybounds[0];

    xp = flint_calloc(xlimbs + 1, sizeof(mp_limb_t));
    yp = flint_calloc(ylimbs + 1, sizeof(mp_limb_t));
    zp = flint_calloc(xlimbs + ylimbs + 1, sizeof(mp_limb_t));

    __nmod_dmpoly_bit_pack(xp, x, xbounds[0], vars, stride + 1, bits);
    __nmod_dmpoly_bit_pack(yp, y, ybounds[0], vars, stride + 1, bits);

    if (xlimbs < ylimbs)
        mpn_mul(zp, yp, ylimbs, xp, xlimbs);
    else
        mpn_mul(zp, xp, xlimbs, yp, ylimbs);

    __nmod_dmpoly_bit_unpack(z, zp, vars, zbounds, stride + 1, bits, mod);

    flint_free(xp);
    flint_free(yp);
    flint_free(zp);
    flint_free(stride);
    flint_free(zbounds);
}


void
_nmod_dmpoly_mul_KS(arr_ptr z, arr_srcptr x, arr_srcptr y, int vars, nmod_t mod)
{
    long xlen, ylen, *xbounds, *ybounds;

    xlen = x->length;
    ylen = y->length;

    if (xlen == 0 || ylen == 0)
    {
        _nmod_dmpoly_zero(z, vars);
        return;
    }

    _nmod_dmpoly_fit_length(z, xlen + ylen - 1, vars);

    if (vars == 1)
    {
        if (xlen >= ylen)
            _nmod_poly_mul_KS(z->coeffs, x->coeffs, xlen,
                                         y->coeffs, ylen, 0, mod);
        else
            _nmod_poly_mul_KS(z->coeffs, y->coeffs, ylen,
                                         x->coeffs, xlen, 0, mod);
    }
    else
    {
        xbounds = flint_calloc(vars, sizeof(long));
        ybounds = flint_calloc(vars, sizeof(long));

        _nmod_dmpoly_rect_hull(xbounds, x, vars);
        _nmod_dmpoly_rect_hull(ybounds, y, vars);

        __nmod_dmpoly_mul_KS(z->coeffs, x->coeffs, y->coeffs,
            xbounds, ybounds, vars, mod);

        flint_free(xbounds);
        flint_free(ybounds);
    }

    z->length = xlen + ylen - 1;
    _nmod_dmpoly_normalise(z, vars);
}

void
nmod_dmpoly_mul_KS(nmod_dmpoly_t z, const nmod_dmpoly_t x,
                                        const nmod_dmpoly_t y)
{
    _nmod_dmpoly_mul_KS(&z->arr, &x->arr, &y->arr, x->vars, x->mod);
}

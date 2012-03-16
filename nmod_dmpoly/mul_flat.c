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

static __inline__ void
__nmod_dmpoly_mul_flat_2D(arr_ptr z, arr_srcptr x, long xlen,
    arr_srcptr y, long ylen, long xybound, long bits, nmod_t mod)
{
    long xi, yi, i, j, zlen, alen, blen, nlimbs;
    mp_limb_t a, b, c, d;
    mp_ptr PZ;
    mp_srcptr PX, PY;
    mp_ptr tmp, t;

    zlen = xlen + ylen - 1;

    if (bits <= FLINT_BITS)
        nlimbs = 1;
    else if (bits <= 2 * FLINT_BITS)
        nlimbs = 2;
    else
        nlimbs = 3;

    tmp = calloc(nlimbs * zlen * xybound, sizeof(mp_limb_t));

    for (i = 0; i < zlen; i++)
        z[i].length = 0;

    for (xi = 0; xi < xlen; xi++)
    {
        alen = x[xi].length;

        for (yi = 0; yi < ylen; yi++)
        {
            blen = y[yi].length;

            PX = (x + xi)->coeffs;
            PY = (y + yi)->coeffs;
            PZ = (z + xi + yi)->coeffs;

            if (nlimbs == 1)
                for (i = 0; i < alen; i++)
                    for (j = 0; j < blen; j++)
                        tmp[xybound*(xi+yi) + i+j] += PX[i] * PY[j];
            else if (nlimbs == 2)
                for (i = 0; i < alen; i++)
                    for (j = 0; j < blen; j++)
                    {
                        t = tmp + 2 * (xybound * (xi+yi) + i+j);
                        umul_ppmm(b, a, PX[i], PY[j]);
                        add_ssaaaa(t[1], t[0], t[1], t[0], b, a);
                    }
            else
                for (i = 0; i < alen; i++)
                    for (j = 0; j < blen; j++)
                    {
                        t = tmp + 3 * (xybound * (xi+yi) + i+j);
                        umul_ppmm(b, a, PX[i], PY[j]);
                        add_sssaaaaaa(t[2], t[1], t[0],
                                      t[2], t[1], t[0], 0, b, a);
                    }

            z[xi+yi].length = FLINT_MAX(z[xi+yi].length, alen + blen - 1);
        }
    }

    for (i = 0; i < zlen; i++)
    {
        _nmod_dmpoly_fit_length(z + i, z[i].length, 1);

        if (nlimbs == 1)
            _nmod_vec_reduce((z + i)->coeffs,
                (mp_srcptr) tmp + xybound * i, z[i].length, mod);
        else if (nlimbs == 2)
            for (j = 0; j < z[i].length; j++)
            {
                a = tmp[2 * (xybound * i + j)];
                b = tmp[2 * (xybound * i + j) + 1];
                c = n_ll_mod_preinv(b, a, mod.n, mod.ninv);
                ((mp_ptr) (z + i)->coeffs)[j] = c;
            }
        else
            for (j = 0; j < z[i].length; j++)
            {
                a = tmp[3 * (xybound * i + j)];
                b = tmp[3 * (xybound * i + j) + 1];
                c = tmp[3 * (xybound * i + j) + 2];
                d = n_lll_mod_preinv(c, b, a, mod.n, mod.ninv);
                ((mp_ptr) (z + i)->coeffs)[j] = d;
            }

        _nmod_dmpoly_normalise(z + i, 1);
    }

    free(tmp);
}

void
__nmod_dmpoly_mul_flat(void * z, const void * x, 
                            const void * y, const long * xbounds,
                const long * ybounds, int vars, nmod_t mod)
{
    long xybound, bits;

    if (vars != 2)
    {
        printf("error: mul_flat only implemented for vars = 2\n");
        abort();
    }

    /* bound inner length */
    xybound = xbounds[1] + ybounds[1] - 1;

    /* bound coefficients */
    bits = FLINT_BIT_COUNT(FLINT_MIN(xbounds[0], ybounds[0]) *
           FLINT_MIN(xbounds[1], ybounds[1])) +
        2*FLINT_BIT_COUNT(mod.n - 1);

    __nmod_dmpoly_mul_flat_2D(z, x, xbounds[0], y, ybounds[0],
        xybound, bits, mod);
}


void
_nmod_dmpoly_mul_flat(arr_ptr z, arr_srcptr x,
                                        arr_srcptr y, int vars, nmod_t mod)
{
    long xlen, ylen, zlen;
    long xbounds[2];
    long ybounds[2];

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
        _nmod_dmpoly_mul_flat(&tmp, x, y, vars, mod);
        _nmod_dmpoly_swap(z, &tmp);
        _nmod_dmpoly_clear(&tmp, vars);
        return;
    }

    xbounds[0] = 0; xbounds[1] = 0;
    _nmod_dmpoly_rect_hull(xbounds, x, vars);

    ybounds[0] = 0; ybounds[1] = 0;
    _nmod_dmpoly_rect_hull(ybounds, y, vars);

    zlen = xlen + ylen - 1;
    _nmod_dmpoly_fit_length(z, zlen, vars);
    if (xlen >= ylen)
        __nmod_dmpoly_mul_flat(z->coeffs, x->coeffs, y->coeffs,
            xbounds, ybounds, vars, mod);
    else
        __nmod_dmpoly_mul_flat(z->coeffs, y->coeffs, x->coeffs,
            ybounds, xbounds, vars, mod);
    z->length = zlen;
    _nmod_dmpoly_normalise(z, vars);
}

void
nmod_dmpoly_mul_flat(nmod_dmpoly_t z, const nmod_dmpoly_t x,
                                    const nmod_dmpoly_t y)
{
    _nmod_dmpoly_mul_flat(&z->arr, &x->arr, &y->arr, x->vars, x->mod);
}

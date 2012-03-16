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

#define MAX_INLINE 20


void
__nmod_dmpoly_mul_2D_classical(arr_ptr z, arr_srcptr x, long xlen,
    arr_srcptr y, long ylen, long xybound, int nlimbs, nmod_t mod)
{
    long xi, yi, i, j, zlen, alen, blen;
    mp_limb_t a, b, c, d;
    mp_ptr PZ;
    mp_srcptr PX, PY;
    mp_ptr tmp, t;

    zlen = xlen + ylen - 1;

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
__nmod_dmpoly_mul_classical(void * z, const void * x, long xlen, const void * y, long ylen, int vars, nmod_t mod)
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

        if (vars == 2 && 0)
        {
            long xmax, ymax, xybound, zlen;

            xmax = 0;
            for (i = 0; i < xlen; i++)
                xmax = FLINT_MAX(xmax, xp[i].length);

            ymax = 0;
            for (i = 0; i < ylen; i++)
                ymax = FLINT_MAX(ymax, yp[i].length);

            xybound = xmax + ymax - 1;
            zlen = xlen + ylen - 1;

            if (xmax < MAX_INLINE && ymax < MAX_INLINE)
            {

                /* Bound number of limbs */


                __nmod_dmpoly_mul_2D_classical(z, x, xlen, y, ylen, xybound, 2, mod);
                return;
            }
        }

        for (i = 0; i < xlen; i++)
            _nmod_dmpoly_mul_classical(zp + i, xp + i, yp, vars - 1, mod);

        if (ylen != 1)
        {
            arr_struct tmp;
            tmp.coeffs = NULL;
            tmp.alloc = 0;
            tmp.length = 0;

            for (i = 0; i < ylen - 1; i++)
                _nmod_dmpoly_mul_classical(zp + xlen + i,
                    yp + i + 1, xp + xlen - 1, vars - 1, mod);

            for (i = 0; i < xlen - 1; i++)
            {
                for (j = 0; j < ylen - 1; j++)
                {
                    _nmod_dmpoly_mul_classical(&tmp, yp + 1 + j, xp + i, vars - 1, mod);
                    _nmod_dmpoly_add(zp + i + j + 1, zp + i + j + 1, &tmp, vars - 1, mod);
                }
            }

            _nmod_dmpoly_clear(&tmp, vars - 1);
        }
    }
}


void
_nmod_dmpoly_mul_classical(arr_ptr z, arr_srcptr x,
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
        _nmod_dmpoly_mul_classical(&tmp, x, y, vars, mod);
        _nmod_dmpoly_swap(z, &tmp);
        _nmod_dmpoly_clear(&tmp, vars);
        return;
    }

    zlen = xlen + ylen - 1;
    _nmod_dmpoly_fit_length(z, zlen, vars);
    if (xlen >= ylen)
        __nmod_dmpoly_mul_classical(z->coeffs, x->coeffs, xlen,
                                            y->coeffs, ylen, vars, mod);
    else
        __nmod_dmpoly_mul_classical(z->coeffs, y->coeffs, ylen,
                                            x->coeffs, xlen, vars, mod);
    z->length = zlen;
    _nmod_dmpoly_normalise(z, vars);
}

void
nmod_dmpoly_mul_classical(nmod_dmpoly_t z, const nmod_dmpoly_t x,
                                    const nmod_dmpoly_t y)
{
    _nmod_dmpoly_mul_classical(&z->arr, &x->arr, &y->arr, x->vars, x->mod);
}

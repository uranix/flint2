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
__nmod_dmpoly_add(void * z, const void * x, long xlen,
                            const void * y, long ylen, int vars, nmod_t mod)
{
    if (vars == 1)
    {
        _nmod_vec_add((mp_ptr) z, (mp_srcptr) x, (mp_srcptr) y, ylen, mod);
        _nmod_vec_set(((mp_ptr) z) + ylen, ((mp_srcptr) x) + ylen, xlen - ylen);
    }
    else
    {
        long i;

        for (i = 0; i < ylen; i++)
            _nmod_dmpoly_add((arr_ptr) z + i, (arr_srcptr) x + i,
                             (arr_srcptr) y + i, vars - 1, mod);

        for (i = ylen; i < xlen; i++)
            _nmod_dmpoly_set((arr_ptr) z + i, (arr_srcptr) x + i, vars - 1);
    }
}

void
_nmod_dmpoly_add(arr_ptr z, arr_srcptr x, arr_srcptr y, int vars, nmod_t mod)
{
    long xlen, ylen;

    xlen = x->length;
    ylen = y->length;

    if (xlen < ylen)
    {
        _nmod_dmpoly_add(z, y, x, vars, mod);
        return;
    }

    if (xlen == 0)
    {
        _nmod_dmpoly_zero(z, vars);
        return;
    }

    _nmod_dmpoly_fit_length(z, xlen, vars);
    __nmod_dmpoly_add(z->coeffs, x->coeffs, xlen, y->coeffs, ylen, vars, mod);
    z->length = xlen;
    _nmod_dmpoly_normalise(z, vars);
}

void
nmod_dmpoly_add(nmod_dmpoly_t z, const nmod_dmpoly_t x, const nmod_dmpoly_t y)
{
    _nmod_dmpoly_add(&z->arr, &x->arr, &y->arr, x->vars, x->mod);
}

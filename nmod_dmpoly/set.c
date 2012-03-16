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
__nmod_dmpoly_set(void * dest, const void * src, long len, int vars)
{
    if (vars == 1)
    {
        _nmod_vec_set((mp_ptr) dest, (mp_srcptr) src, len);
    }
    else
    {
        long i;
        for (i = 0; i < len; i++)
            _nmod_dmpoly_set((arr_ptr) dest + i,
                            (arr_srcptr) src + i, vars - 1);
    }
}

void
_nmod_dmpoly_set(arr_ptr dest, arr_srcptr src, int vars)
{
    if (src->length == 0)
    {
        _nmod_dmpoly_zero(dest, vars);
    }
    else
    {
        _nmod_dmpoly_fit_length(dest, src->length, vars);
        __nmod_dmpoly_set(dest->coeffs, src->coeffs, src->length, vars);
        dest->length = src->length;
    }
}

void
nmod_dmpoly_set(nmod_dmpoly_t dest, const nmod_dmpoly_t src)
{
    _nmod_dmpoly_set(&dest->arr, &src->arr, src->vars);
}

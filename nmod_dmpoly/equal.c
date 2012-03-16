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

int
__nmod_dmpoly_equal(const void * p1, const void * p2, long len, int vars)
{
    long i;

    if (vars == 1)
        return _nmod_vec_equal((mp_ptr) p1, (mp_srcptr) p2, len);

    for (i = 0; i < len; i++)
        if (!_nmod_dmpoly_equal((arr_srcptr) p1 + i,
                                (arr_srcptr) p2 + i, vars - 1))
            return 0;

    return 1;
}

int
_nmod_dmpoly_equal(arr_srcptr p1, arr_srcptr p2, int vars)
{
    if (p1->length != p2->length)
        return 0;

    return __nmod_dmpoly_equal(p1->coeffs, p2->coeffs, p2->length, vars);
}

int
nmod_dmpoly_equal(const nmod_dmpoly_t p1, const nmod_dmpoly_t p2)
{
    return _nmod_dmpoly_equal(&p1->arr, &p2->arr, p2->vars);
}

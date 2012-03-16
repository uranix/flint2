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
_nmod_dmpoly_randtest(arr_ptr poly, flint_rand_t state, long maxlen,
    int vars, nmod_t mod)
{
    long i, n = n_randint(state, maxlen + 1);

    _nmod_dmpoly_fit_length(poly, n, vars);

    if (vars == 1)
    {
        mp_ptr ptr = (mp_ptr) poly->coeffs;
        _nmod_vec_randtest(ptr, state, n, mod);
    }
    else
    {
        arr_ptr ptr = (arr_ptr) poly->coeffs;
        for (i = 0; i < n; i++)
            _nmod_dmpoly_randtest(ptr + i, state, maxlen, vars - 1, mod);
    }

    poly->length = n;
    _nmod_dmpoly_normalise(poly, vars);
}

void
nmod_dmpoly_randtest(nmod_dmpoly_t poly, flint_rand_t state, long maxlen)
{
    _nmod_dmpoly_randtest(&poly->arr, state, maxlen, poly->vars, poly->mod);
}

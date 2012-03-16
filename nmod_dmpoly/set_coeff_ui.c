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
_nmod_dmpoly_set_coeff_ui(arr_ptr poly, long * pos, mp_limb_t c, int vars)
{
    long i, n = pos[0];

    _nmod_dmpoly_fit_length(poly, n + 1, vars);

    if (vars == 1)
    {
        mp_ptr ptr = (mp_ptr) poly->coeffs;
        ptr[n] = c;
        for (i = poly->length; i < n; i++)
            ptr[i] = 0;
    }
    else
    {
        arr_ptr ptr = (arr_ptr) poly->coeffs;
        _nmod_dmpoly_set_coeff_ui(ptr + n, pos + 1, c, vars - 1);
        for (i = poly->length; i < n; i++)
            _nmod_dmpoly_zero(ptr + i, vars - 1);
    }

    poly->length = FLINT_MAX(poly->length, n + 1);
    _nmod_dmpoly_normalise(poly, vars);
}

void
nmod_dmpoly_set_coeff_ui(nmod_dmpoly_t poly, long * pos, mp_limb_t c)
{
    c = n_mod2_preinv(c, poly->mod.n, poly->mod.ninv);
    _nmod_dmpoly_set_coeff_ui(&poly->arr, pos, c, poly->vars);
}

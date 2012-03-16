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

mp_limb_t
__nmod_dmpoly_evaluate_nmod(const void * poly, long len, mp_srcptr x,
    int vars, nmod_t mod)
{
    if (vars == 1)
    {
        return _nmod_poly_evaluate_nmod((mp_srcptr) poly, len, x[0], mod);
    }
    else
    {
        long i = len - 1;
        mp_limb_t t, val;

        val = _nmod_dmpoly_evaluate_nmod((arr_srcptr) poly + i,
                    x + 1, vars - 1, mod);
        i--;

        for ( ; i >= 0; i--)
        {
            val = n_mulmod2_preinv(val, x[0], mod.n, mod.ninv);
            t = _nmod_dmpoly_evaluate_nmod((arr_srcptr) poly + i,
                    x + 1, vars - 1, mod);
            val = n_addmod(val, t, mod.n);
        }

        return val;
    }
}

mp_limb_t
_nmod_dmpoly_evaluate_nmod(arr_srcptr poly, mp_srcptr x, int vars, nmod_t mod)
{
    if (poly->length == 0)
        return 0;

    return __nmod_dmpoly_evaluate_nmod(poly->coeffs,
        poly->length, x, vars, mod);
}

mp_limb_t
nmod_dmpoly_evaluate_nmod(const nmod_dmpoly_t poly, mp_srcptr x)
{
    return _nmod_dmpoly_evaluate_nmod(&poly->arr, x, poly->vars, poly->mod);
}
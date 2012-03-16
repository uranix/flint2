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
_nmod_dmpoly_get_coeff_ui(arr_srcptr poly, long * pos, int vars)
{
    if (pos[0] >= poly->length)
        return 0UL;

    if (vars == 1)
    {
        return ((mp_srcptr) poly->coeffs)[pos[0]];
    }
    else
    {
        return _nmod_dmpoly_get_coeff_ui((arr_srcptr) poly->coeffs + pos[0],
            pos + 1, vars - 1);
    }
}

mp_limb_t
nmod_dmpoly_get_coeff_ui(const nmod_dmpoly_t poly, long * pos)
{
    return _nmod_dmpoly_get_coeff_ui(&poly->arr, pos, poly->vars);
}

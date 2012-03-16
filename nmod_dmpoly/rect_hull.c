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

/* Computes the maximum length in each variable, assuming that
   the lengths vector has been pre-initialised to zero.
   Returns the total length of all polynomials at the bottom level. */

long
__nmod_dmpoly_rect_hull(long * bounds, const void * poly, long len, int vars)
{
    bounds[0] = FLINT_MAX(bounds[0], len);

    if (vars == 1)
    {
        return len;
    }
    else if (vars == 2)  /* inline */
    {
        long total, l, i;

        for (i = total = 0; i < len; i++)
        {
            l = ((arr_srcptr) poly)[i].length;
            bounds[1] = FLINT_MAX(bounds[1], l);
            total += l;
        }

        return total;
    }
    else
    {
        long total, i;

        for (i = total = 0; i < len; i++)
            total += _nmod_dmpoly_rect_hull(bounds + 1,
                    (arr_srcptr) poly + i, vars - 1);

        return total;
    }
}

long
_nmod_dmpoly_rect_hull(long * bounds, arr_srcptr poly, int vars)
{
    return __nmod_dmpoly_rect_hull(bounds, poly->coeffs, poly->length, vars);
}

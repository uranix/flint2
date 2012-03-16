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

void _nmod_dmpoly_fit_length(arr_ptr poly, long len, int vars)
{
    if (len > poly->alloc)
    {
        long i, alloc, new_alloc;

        alloc = poly->alloc;
        new_alloc = FLINT_MAX(len, 2 * alloc);
        poly->alloc = new_alloc;

        if (vars == 1)
        {
            poly->coeffs = flint_realloc(poly->coeffs,
                new_alloc * sizeof(mp_limb_t));

            for (i = alloc; i < new_alloc; i++)
                ((mp_ptr) poly->coeffs)[i] = 0;
        }
        else
        {
            arr_ptr ptr;

            poly->coeffs = ptr = flint_realloc(poly->coeffs,
                new_alloc * sizeof(arr_struct));

            for (i = alloc; i < new_alloc; i++)
            {
                ptr[i].coeffs = NULL;
                ptr[i].alloc = 0;
                ptr[i].length = 0;
            }
        }
    }
}

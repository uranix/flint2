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
_nmod_dmpoly_clear(arr_struct * poly, int vars)
{
    if (poly->alloc != 0)
    {
        if (vars > 1)
        {
            long i, alloc = poly->alloc;
            arr_ptr ptr = (arr_ptr) poly->coeffs;

            for (i = 0; i < alloc; i++)
                _nmod_dmpoly_clear(ptr + i, vars - 1);
        }

        flint_free(poly->coeffs);
    }
}

void
nmod_dmpoly_clear(nmod_dmpoly_t poly)
{
    _nmod_dmpoly_clear(&poly->arr, poly->vars);
}

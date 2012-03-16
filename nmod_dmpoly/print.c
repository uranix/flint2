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
__nmod_dmpoly_print(const void * poly, long len, int vars)
{
    long i;

    printf("[");

    if (vars == 1)
    {
        mp_srcptr ptr = (mp_srcptr) poly;

        for (i = 0; i < len; i++)
        {
            printf("%lu", ptr[i]);
            if (i + 1 < len)
                printf(", ");
        }
    }
    else
    {
        arr_srcptr ptr = (arr_srcptr) poly;

        for (i = 0; i < len; i++)
        {
            _nmod_dmpoly_print(ptr + i, vars - 1);
            if (i + 1 < len)
                printf(", ");
        }
    }

    printf("]");
}

void
_nmod_dmpoly_print(arr_srcptr poly, int vars)
{
    __nmod_dmpoly_print(poly->coeffs, poly->length, vars);
}

void
nmod_dmpoly_print(const nmod_dmpoly_t poly)
{
    _nmod_dmpoly_print(&poly->arr, poly->vars);
}

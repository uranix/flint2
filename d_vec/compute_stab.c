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

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "flint.h"
#include "d_vec.h"

#define PI 3.141592653589793
/*
   computes nh/2 values in stab where nh = n/2, namely
   sin(0*pi/nh), sin(1*pi/nh), ..., sin((nh-1)*pi/nh
*/
void
_d_vec_compute_stab(double * stab, long n)
{
    long nh, i;
    
    while (n >= 1)
    {
        stab[0] = 0.0;

        nh = n/2;

        for (i = 1; i < nh/2; i++)
           stab[i] = sin((PI*i)/nh);

        stab += nh/2;
        n /= 8;
    }
}

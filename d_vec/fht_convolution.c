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
#include "flint.h"
#include "d_vec.h"

/*
   requires nh/2 values in stab where nh = n/2, namely
   sin(0*pi/nh), sin(1*pi/nh), ..., sin((nh-1)*pi/nh
*/
void
_d_vec_fht_convolution(double * in1, double * in2, double * stab, long n)
{
    long i, nh;
    double t;
    double * p1, * p2;

    _d_vec_fht_dif_rec(in2, stab, 1, n);
    _d_vec_fht_dif_rec(in1, stab, 1, n);

    nh = n/2;

    in1[0] *= in2[0];
    if (n > 1)
    {
        in1[1] *= in2[1];
        n = 2;
    } else n = 1;

    p1 = in1 + 2;
    p2 = in2 + 2;
    while (n <= nh)
    {
        for (i = 0; i < n/2; i++)
        {
           t = (p2[i]*(p1[i] + p1[n - i - 1]) + p2[n - i - 1]*(p1[i] - p1[n - i - 1]))*0.5;
           p1[n - i - 1] = (p2[n - i - 1]*(p1[n - i - 1] + p1[i]) + p2[i]*(p1[n - i - 1] - p1[i]))*0.5;
           p1[i] = t;
        }
        p1 += n;
        p2 += n;
        n *= 2;
    }
    
    _d_vec_fht_dit_rec(in1, stab, 1, n);
    t = 1.0/n;
    for (i = 0; i < n; i++)
        in1[i] *= t;
}

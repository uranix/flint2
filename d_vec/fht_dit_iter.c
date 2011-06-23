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
   they are expected to be separated by step in the array stab
*/
void
_d_vec_fht_dit_iter(double * in, double * stab, long step, long n)
{
    long nh, i, send, n2 = 2;
    double t;
    double * s1, * end = in + n;

    step *= (n/2);

    while (n2 <= n)
    {
        nh = n2/2;
        s1 = in;
        send = (step*nh)/2;

        while (s1 < end)
        {
            for (i = 1; i < nh/2; i++)
            {
                t = s1[i + nh]*stab[send - i*step] + s1[n2 - i] *stab[i*step];
                s1[n2 - i] = s1[i + nh]*stab[i*step] - s1[n2 - i]*stab[send - i*step];
                s1[i + nh] = t;
            }
            
            for (i = 0; i < nh; i++)
            {
                t = s1[i] - s1[i + nh];
                s1[i] += s1[i + nh];
                s1[i + nh] = t;
            }

            s1 += n2;
        }
        n2 *= 2;
        step /= 2;
    }
}

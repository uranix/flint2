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
_d_vec_fht_dit_rec(double * in, double * stab, long step, long n)
{
    long nh, i, send;
    double t;

    if (n <= D_VEC_FHT_REC_CUTOFF)
    {
        _d_vec_fht_dit_iter(in, stab, step, n);
        return;
    }

    nh = n/2;
    send = (step*nh)/2;

    if (step == 4)
    {
        _d_vec_fht_dit_rec(in     , stab + send, 1, nh);
        _d_vec_fht_dit_rec(in + nh, stab + send, 1, nh);
    } 
    else 
    {
        _d_vec_fht_dit_rec(in     , stab, 2*step, nh);
        _d_vec_fht_dit_rec(in + nh, stab, 2*step, nh);

    }

    for (i = 1; i < nh/2; i++)
    {
        t = in[i + nh]*stab[send - i*step] + in[n - i] *stab[i*step];
        in[n - i] = in[i + nh]*stab[i*step] - in[n - i]*stab[send - i*step];
        in[i + nh] = t;
    }

    for (i = 0; i < nh; i++)
    {
        t = in[i] - in[i + nh];
        in[i] += in[i + nh];
        in[i + nh] = t;
    }
}

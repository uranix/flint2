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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "nmod_dmpoly.h"


static __inline__ void
__nmod_dmpoly_pow_binexp(arr_ptr res, arr_srcptr poly, long len,
        ulong e, int vars, nmod_t mod)
{
    ulong bit = ~((~0UL) >> 1);
    long rlen;
    long alloc = (long) e * (len - 1) + 1;
    arr_struct tmp;
    arr_ptr v, R, S, T;

    _nmod_dmpoly_init(&tmp, vars);
    _nmod_dmpoly_fit_length(&tmp, alloc, vars);
    v = tmp.coeffs;

    while ((bit & e) == 0UL)
        bit >>= 1;
    bit >>= 1;

    /*
       Trial run without any polynomial arithmetic to determine the parity 
       of the number of swaps;  then set R and S accordingly
     */
    {
        unsigned int swaps = 0U;
        ulong bit2 = bit;
        if ((bit2 & e))
            swaps = ~swaps;
        while (bit2 >>= 1)
            if ((bit2 & e) == 0UL)
                swaps = ~swaps;
        
        if (swaps == 0U)
        {
            R = res;
            S = v;
        }
        else
        {
            R = v;
            S = res;
        }
    }

    /*
       We unroll the first step of the loop, referring to {poly, len}
     */
    
    __nmod_dmpoly_mul(R, poly, len, poly, len, vars, mod);
    rlen = 2 * len - 1;

    if ((bit & e))
    {
        __nmod_dmpoly_mul(S, R, rlen, poly, len, vars, mod);
        rlen += len - 1;
        T = R;
        R = S;
        S = T;
    }
    
    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            __nmod_dmpoly_mul(S, R, rlen, R, rlen, vars, mod);
            rlen += rlen - 1;
            __nmod_dmpoly_mul(R, S, rlen, poly, len, vars, mod);
            rlen += len - 1;
        }
        else
        {
            __nmod_dmpoly_mul(S, R, rlen, R, rlen, vars, mod);
            rlen += rlen - 1;
            T = R;
            R = S;
            S = T;
        }
    }

    _nmod_dmpoly_clear(&tmp, vars);
}


void
__nmod_dmpoly_pow(void * z, const void * x, long xlen,
        ulong exp, int vars, nmod_t mod)
{
    if (vars == 1)
        _nmod_poly_pow((mp_ptr) z, (mp_srcptr) x, xlen, exp, mod);
    else
        __nmod_dmpoly_pow_binexp((arr_ptr) z,
            (arr_srcptr) x, xlen, exp, vars, mod);
}


void
_nmod_dmpoly_pow(arr_ptr z, arr_srcptr x, ulong exp, int vars, nmod_t mod)
{
    long xlen, zlen;

    xlen = x->length;

    if (exp == 0)
    {
        _nmod_dmpoly_set_ui(z, 1UL, vars);
        return;
    }

    if (xlen == 0)
    {
        _nmod_dmpoly_zero(z, vars);
        return;
    }

    if (exp == 1)
    {
        _nmod_dmpoly_set(z, x, vars);
        return;
    }

    if (z == x)
    {
        arr_struct tmp;
        _nmod_dmpoly_init(&tmp, vars);
        _nmod_dmpoly_pow(&tmp, x, exp, vars, mod);
        _nmod_dmpoly_swap(z, &tmp);
        _nmod_dmpoly_clear(&tmp, vars);
        return;
    }

    zlen = (long) exp * (xlen - 1) + 1;
    _nmod_dmpoly_fit_length(z, zlen, vars);
    __nmod_dmpoly_pow(z->coeffs, x->coeffs, xlen, exp, vars, mod);
    z->length = zlen;
    _nmod_dmpoly_normalise(z, vars);
}

void
nmod_dmpoly_pow(nmod_dmpoly_t z, const nmod_dmpoly_t x, ulong exp)
{
    _nmod_dmpoly_pow(&z->arr, &x->arr, exp, x->vars, x->mod);
}

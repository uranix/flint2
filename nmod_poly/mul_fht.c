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

    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "d_vec.h"

void
_nmod_poly_mul_fht(mp_ptr out, mp_srcptr in1, long len1,
                  mp_srcptr in2, long len2, mp_bitcnt_t bits, nmod_t mod)
{
    long len_out = len1 + len2 - 1, n, i;
    mp_bitcnt_t bits1, bits2, loglen;
    double * d1, * d2, * stab;

    if (bits == 0)
    {
        bits1  = _nmod_vec_max_bits(in1, len1);
        bits2  = (in1 == in2) ? bits1 : _nmod_vec_max_bits(in2, len2);
        loglen = FLINT_BIT_COUNT(len2);
        
        bits = bits1 + bits2 + loglen;
    }

    loglen = FLINT_BIT_COUNT(len_out);
    n = (1L<<loglen);
    
    d1 = _d_vec_init(n);
    d2 = _d_vec_init(n);
    stab = _d_vec_init((5*n)/16);
    
    _d_vec_compute_stab(stab, n);
      
    for (i = 0; i < len1; i++)
        d1[i] = (double) in1[i];    
    for ( ; i < n; i++)
        d1[i] = 0.0;

    for (i = 0; i < len2; i++)
        d2[i] = (double) in2[i];    
    for ( ; i < n; i++)
        d2[i] = 0.0;

    _d_vec_fht_convolution(d1, d2, stab, n);

    _d_vec_free(d2);
    _d_vec_free(stab);

    for (i = 0; i < len_out; i++)
        NMOD_RED(out[i], (ulong) (d1[i] + 0.5), mod);

    _d_vec_free(d1);
}

void
nmod_poly_mul_fht(nmod_poly_t res,
                 const nmod_poly_t poly1, const nmod_poly_t poly2,
                 mp_bitcnt_t bits)
{
    long len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_fht(temp->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, bits,
                              poly1->mod);
        else
            _nmod_poly_mul_fht(temp->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length, bits,
                              poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_fht(res->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, bits,
                              poly1->mod);
        else
            _nmod_poly_mul_fht(res->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length, bits,
                              poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}

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

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_lcm(fmpz * res, const fmpz * poly1, long len1,
                           const fmpz * poly2, long len2)
{
    fmpz *v;
    long lenv;

    v = _fmpz_vec_init(len2);
    _fmpz_poly_gcd(v, poly1, len1, poly2, len2);
    for (lenv = len2 - 1; (lenv >= 0) && !v[lenv]; lenv--) ;
    lenv++;

    if (lenv == 1 && fmpz_is_one(v))
    {
        _fmpz_poly_mul(res, poly1, len1, poly2, len2);
    }
    else
    {
        fmpz *w;
        long lenw = len1 - lenv + 1;

        w = _fmpz_vec_init(lenw);
        _fmpz_poly_divides(w, poly1, len1, v, lenv);
        if (lenw >= len2)
        {
            _fmpz_poly_mul(res, w, lenw, poly2, len2);
        }
        else
        {
            _fmpz_poly_mul(res, poly2, len2, w, lenw);
        }
        _fmpz_vec_zero(res + (lenw + len2 - 1), len1 - lenw);
        _fmpz_vec_clear(w, lenw);
    }
    if (fmpz_sgn(res + (len1 + len2 - lenv - 1)) < 0)
    {
        _fmpz_vec_neg(res, res, len1 + len2 - lenv);
    }
    _fmpz_vec_clear(v, len2);
}

void
fmpz_poly_lcm(fmpz_poly_t res, const fmpz_poly_t poly1,
                               const fmpz_poly_t poly2)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;

        fmpz_poly_init(t);
        fmpz_poly_lcm(t, poly1, poly2);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    fmpz_poly_fit_length(res, len1 + len2 - 1);
    if (len1 >= len2)
    {
        _fmpz_poly_lcm(res->coeffs, poly1->coeffs, len1, 
                                    poly2->coeffs, len2);
    }
    else
    {
        _fmpz_poly_lcm(res->coeffs, poly2->coeffs, len2, 
                                    poly1->coeffs, len1);
    }
    _fmpz_poly_set_length(res, len1 + len2 - 1);
    _fmpz_poly_normalise(res);
    
    return;
}


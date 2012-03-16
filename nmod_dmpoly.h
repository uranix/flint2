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

#ifndef NMOD_DMPOLY_H
#define NMOD_DMPOLY_H

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

typedef struct
{
    void * coeffs;
    long length;
    long alloc;
}
arr_struct;

typedef arr_struct * arr_ptr;
typedef const arr_struct * arr_srcptr;

typedef struct
{
    arr_struct arr;
    nmod_t mod;
    int vars;
} nmod_dmpoly_struct;

typedef nmod_dmpoly_struct nmod_dmpoly_t[1];

/* Memory management *********************************************************/

static __inline__ void _nmod_dmpoly_init(arr_ptr poly, int vars)
{
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void nmod_dmpoly_init(nmod_dmpoly_t poly, int vars, mp_limb_t mod);

void _nmod_dmpoly_clear(arr_struct * poly, int vars);

void nmod_dmpoly_clear(nmod_dmpoly_t poly);

void _nmod_dmpoly_fit_length(arr_ptr poly, long len, int vars);

void _nmod_dmpoly_normalise(arr_ptr poly, int vars);

/* Basic manipulation ********************************************************/

void __nmod_dmpoly_set(void * dest, const void * src, long len, int vars);

void _nmod_dmpoly_set(arr_ptr dest, arr_srcptr src, int vars);

void nmod_dmpoly_set(nmod_dmpoly_t dest, const nmod_dmpoly_t src);

static __inline__ void _nmod_dmpoly_swap(arr_ptr f, arr_ptr g)
{
    arr_struct tmp = *f;
    *f = *g;
    *g = tmp;
}

void nmod_dmpoly_swap(nmod_dmpoly_t a, nmod_dmpoly_t b);

void __nmod_dmpoly_zero(void * poly, long len, int vars);

void _nmod_dmpoly_zero(arr_struct * poly, int vars);

void nmod_dmpoly_zero(nmod_dmpoly_t poly);

void _nmod_dmpoly_set_ui(arr_ptr poly, mp_limb_t c, int vars);

void nmod_dmpoly_set_ui(nmod_dmpoly_t poly, mp_limb_t c);

/* IO ************************************************************************/

void __nmod_dmpoly_print(const void * poly, long len, int vars);

void _nmod_dmpoly_print(arr_srcptr poly, int vars);

void nmod_dmpoly_print(const nmod_dmpoly_t poly);

/* Random generation *********************************************************/

void _nmod_dmpoly_randtest(arr_ptr poly, flint_rand_t state, long maxlen,
    int vars, nmod_t mod);

void nmod_dmpoly_randtest(nmod_dmpoly_t poly, flint_rand_t state, long maxlen);

/* Comparison ****************************************************************/

int __nmod_dmpoly_equal(const void * p1, const void * p2, long len, int vars);

int _nmod_dmpoly_equal(arr_srcptr p1, arr_srcptr p2, int vars);

int nmod_dmpoly_equal(const nmod_dmpoly_t p1, const nmod_dmpoly_t p2);

/* Coefficient manipulation **************************************************/

mp_limb_t _nmod_dmpoly_get_coeff_ui(arr_srcptr poly, long * pos, int vars);

mp_limb_t nmod_dmpoly_get_coeff_ui(const nmod_dmpoly_t poly, long * pos);

void _nmod_dmpoly_set_coeff_ui(arr_ptr poly, long * pos, mp_limb_t c, int vars);

void nmod_dmpoly_set_coeff_ui(nmod_dmpoly_t poly, long * pos, mp_limb_t c);

/* Evaluation ****************************************************************/

mp_limb_t __nmod_dmpoly_evaluate_nmod(const void * poly, long len, mp_srcptr x,
    int vars, nmod_t mod);

mp_limb_t
_nmod_dmpoly_evaluate_nmod(arr_srcptr poly, mp_srcptr x, int vars, nmod_t mod);

mp_limb_t nmod_dmpoly_evaluate_nmod(const nmod_dmpoly_t poly, mp_srcptr x);

/* Bounds ********************************************************************/

long __nmod_dmpoly_rect_hull(long * bounds, const void * poly, long len, int vars);

long _nmod_dmpoly_rect_hull(long * bounds, arr_srcptr poly, int vars);

/* Arithmetic ****************************************************************/

void __nmod_dmpoly_add(void * z, const void * x, long xlen,
                            const void * y, long ylen, int vars, nmod_t mod);
void _nmod_dmpoly_add(arr_ptr z, arr_srcptr x, arr_srcptr y, int vars, nmod_t mod);
void nmod_dmpoly_add(nmod_dmpoly_t z, const nmod_dmpoly_t x, const nmod_dmpoly_t y);

void __nmod_dmpoly_mul(void * z, const void * x, long xlen,
                            const void * y, long ylen, int vars, nmod_t mod);
void _nmod_dmpoly_mul(arr_ptr z, arr_srcptr x,
                                        arr_srcptr y, int vars, nmod_t mod);
void nmod_dmpoly_mul(nmod_dmpoly_t z, const nmod_dmpoly_t x,
                                    const nmod_dmpoly_t y);

void __nmod_dmpoly_mul_flat(void * z, const void * x, 
                            const void * y, const long * xbounds,
                const long * ybounds, int vars, nmod_t mod);
void _nmod_dmpoly_mul_flat(arr_ptr z, arr_srcptr x,
                                        arr_srcptr y, int vars, nmod_t mod);
void nmod_dmpoly_mul_flat(nmod_dmpoly_t z, const nmod_dmpoly_t x,
                                    const nmod_dmpoly_t y);

void __nmod_dmpoly_mul_recursive(void * z, const void * x, long xlen,
                        const void * y, long ylen, int vars, nmod_t mod);
void _nmod_dmpoly_mul_recursive(arr_ptr z, arr_srcptr x,
                                        arr_srcptr y, int vars, nmod_t mod);
void nmod_dmpoly_mul_recursive(nmod_dmpoly_t z, const nmod_dmpoly_t x,
        const nmod_dmpoly_t y);

void __nmod_dmpoly_mul_classical(void * z, const void * x, long xlen,
        const void * y, long ylen, int vars, nmod_t mod);
void _nmod_dmpoly_mul_classical(arr_ptr z, arr_srcptr x,
                                        arr_srcptr y, int vars, nmod_t mod);
void nmod_dmpoly_mul_classical(nmod_dmpoly_t z, const nmod_dmpoly_t x,
        const nmod_dmpoly_t y);

void __nmod_dmpoly_mul_KS(arr_ptr z, arr_srcptr x, arr_srcptr y,
    const long * xbounds, const long * ybounds, int vars, nmod_t mod);
void _nmod_dmpoly_mul_KS(arr_ptr z, arr_srcptr x, arr_srcptr y, int vars, nmod_t mod);
void nmod_dmpoly_mul_KS(nmod_dmpoly_t z, const nmod_dmpoly_t x, const nmod_dmpoly_t y);

void __nmod_dmpoly_pow(void * z, const void * x, long xlen, ulong exp, int vars, nmod_t mod);
void _nmod_dmpoly_pow(arr_ptr z, arr_srcptr x, ulong exp, int vars, nmod_t mod);
void nmod_dmpoly_pow(nmod_dmpoly_t z, const nmod_dmpoly_t x, ulong exp);

#endif

/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

typedef struct cf_ball {
    fmpz_t a, b, e, d
} fmpq_ball_struct;

typedef fmpq_ball fmpq_ball_struct[1];


typedef fmpz_mat22 {
    fmpz_t _11, _12, _21, _22;
    int det;    /* 0,1,or,-1: 0 -> dont know, 1 -> 1, -1 -> -1 */
} fmpz_ball_struct;

typedef fmpz_mat22_t fmpz_mat22_struct[1];


void fmpz_mat22_init(fmpz_mat22_t M)
{
    fmpz_init(M._11);
    fmpz_init(M._12);
    fmpz_init(M._21);
    fmpz_init(M._22);
    M.det = 0;
}

void fmpz_mat22_clear(fmpz_mat22_t M)
{
    fmpz_clear(M._11);
    fmpz_clear(M._12);
    fmpz_clear(M._21);
    fmpz_clear(M._22);
}

/* M = M.N */
void fmpz_mat22_rmul(fmpz_mat22_t M, const fmpz_mat22_t N)
{
    fmpz_t a, b, c, d;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(e);
    fmpz_init(d);
    fmpz_mul(a, M._11, N._11); fmpz_addmul(a, M._12, N._21);
    fmpz_mul(b, M._11, N._12); fmpz_addmul(b, M._12, N._22);
    fmpz_mul(c, M._21, N._11); fmpz_addmul(c, M._22, N._21);
    fmpz_mul(d, M._21, N._12); fmpz_addmul(d, M._22, N._22);
    fmpz_swap(M._11, a.data);
    fmpz_swap(M._12, b.data);
    fmpz_swap(M._21, c.data);
    fmpz_swap(M._22, d.data);
    M->det *= N->det;
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);
}

/* M = M.[q 1; 1 0] */
void fmpz_mat22_rmul_elem(fmpz_mat22_t M, const fmpz_t q)
{
    fmpz_addmul(M._21, M._11, q);
    fmpz_addmul(M._22, M._21, q);
    fmpz_swap(M._11, M._12);
    fmpz_swap(M._21, M._22);
    M->det *= -1;
}



static fmpq_ball_init(fmpq_ball_t x)
{
    fmpz_init(x.a);
    fmpz_init(x.b);
    fmpz_init(x.e);
    fmpz_init(x.d);
}

static fmpq_ball_clear(fmpq_ball_t x)
{
    fmpz_clear(x.a);
    fmpz_clear(x.b);
    fmpz_clear(x.e);
    fmpz_clear(x.d);
}

/* x is canonical and bounded away from 1, i.e. 1 < x */
static fmpq_ball_is_ok(const fmpq_ball_t x)
{
    fmpz_t temp;
    int r;

    if (fmpz_sgn(x.a) <= 0)
        return false;
    if (fmpz_sgn(x.b) <= 0)
        return false;
    if (fmpz_sgn(x.e) < 0)
        return false;
    if (fmpz_sgn(x.d) < 0)
        return false;

    if (fmpz_cmp(x.b, x.a) >= 0)
        return false;

    fmpz_init(temp);
    fmpz_add(temp, x.b, x.e);
    r = fmpz_cmp(temp, temp, x.a) < 0;
    fmpz_clear(temp);
    return r;
}

/* y = m^-1(x) */
static fmpq_ball_apply_mat22_inv(
    fmpq_ball_t y
    const fmpq_mat_t m,
    const fmpq_ball_t x)
{
    if (m.det == 1)
    {
        fmpz_mul(y.a, x.a, m._22);
        fmpz_submul(y.a, x.b, m._12);
        fmpz_submul(y.a, x.e, m._12);
        fmpz_mul(y.b, x.b, m._11);
        fmpz_submul(y.b, x.a, m._21);
        fmpz_submul(y.b, x.d, m._21);
    }
    else
    {
        FLINT_ASSERT(m.det == -1);
        fmpz_mul(y.a, x.b, m._12);
        fmpz_submul(y.a, x.a, m._22);
        fmpz_submul(y.a, x.d, m._22);
        fmpz_mul(y.b, x.a, m._21);
        fmpz_submul(y.b, x.b, m._11);
        fmpz_submul(y.b, x.e, m._11);
    }
    fmpz_mul(y.e, x.e, m._11);
    fmpz_addmul(y.e, x.d, m._21);
    fmpz_mul(y.d, x.e, m._12);
    fmpz_addmul(y.d, x.d, m._22);
}

/* y = [q 1; 1 0]^-1(x) */
static fmpq_ball_apply_mat22_inv_elem(
    fmpq_ball_t y
    const fmpz_t q,
    const fmpq_ball_t x)
{
    fmpz_set(y.a, x.b);
    fmpz_set(y.b, x.a);
    fmpz_submul(y.b, x.b, q);
    fmpz_submul(y.b, x.e, q);
    fmpz_mul(y.e, x.e, q);
    fmpz_add(y.e, x.d);
    fmpz_set(y.d, x.e);
}


/* the interface is entirely mutable and the terms are streamed to s */
fmpq_ball_get_cfrac(
    fmpz_poly_t s,
    fmpz_mat22_t M,
    fmpq_ball_t x,
    const slong limit)
{
    flint_bitcnt_t k;
    fmpq_ball_t y;
    fmpz_mat22_t N;

    FLINT_ASSERT(limit >= 0);

    fmpq_ball_init(y);
    fmpz_mat22_init(M)

again:
    FLINT_ASSERT(count <= limit);
    FLINT_ASSERT(fmpq_ball_is_ok(x));

    if (s->length >= limit)       
        goto cleanup;

    k = fmpz_bits(a);

    if k > FLINT_BITS
        goto split;

gauss:

    fmpz_fdiv_qr(q, r, x.a, x.b);
    fmpq_ball_apply_mat22_inv_elem(y, q, x);
    if (!fmpq_ball_is_ok(y))
        goto cleanup;
    fmpq_ball_swap(x, y);
    fmpz_mat22_rmul_elem(M, q);
    fmpz_poly_fit_length(s, s->length + 1);
    fmpz_swap(s->coeffs + s->length, q);
    s->length++;
    goto again;
    

split:
    k = k/2;

    fmpz_fdiv_q_2exp(y.a, x.a, k);
    fmpz_fdiv_q_2exp(y.b, x.b, k);
    fmpz_fdiv_q_2exp(y.e, x.e, k);
    fmpz_fdiv_q_2exp(y.d, x.d, k);
    fmpz_add_ui(y.e, y.e, 2);
    fmpz_add_ui(y.d, y.d, 2);
    if (!fmpq_ball_is_ok(y))
        goto gauss;
    fmpz_mat22_one(N);
    old_length = s->length;
    fmpq_ball_get_cfrac(s, N, y, limit);
    if (old_length <= s->length)
    {
        FLINT_ASSERT(old_length == s->length);
        FLINT_ASSERT(fmpz_mat22_is_one(N));
        goto gauss;
    }
    fmpq_ball_apply_mat22_inv_elem(x, N, y);
    fmpz_mat22_rmul(M, N);
    goto again;

cleanup:
    FLINT_ASSERT(fmpz_poly_length(s) <= limit);
    fmpq_ball_clear(y);
    fmpz_mat22_clear(N);
    return;
}


slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t f, slong limit)
{
    fmpq_ball_t x;
    fmpz_poly_t s;
    fmpz_mat22_t M;

    FLINT_ASSERT(fmpq_is_canonical(f));

    if (limit <= 0)
        return;

    fmpz_mat22_init(M);
    fmpz_mat22_one(M);

    fmpq_ball_init(x);
    fmpz_abs(x.a, fmpq_numref(f));
    fmpz_abs(x.b, fmpq_denref(f));
    fmpz_zero(x.e);
    fmpz_zero(x.d);

    fmpz_poly_init(s);
    s->length = 0;

    cmp = fmpz_cmp(x.a, x.b);   
    if (cmp > 0)
    {
        fmpq_ball_get_cfrac(s, M, x, limit);
    }
    else if (fmpz
    {
        fmpz_poly_fit_length(s, 1);
        fmpz_zero(s->coeffs + 0);
        s->length = 1;
        fmpz_swap(x.a, x.b);
        fmpq_ball_get_cfrac(s, M, x, limit);
    }

    FLINT_ASSERT(fmpz_is_zero(x.e));
    FLINT_ASSERT(fmpz_is_zero(x.d));
    FLINT_ASSERT(fmpq_ball_is_ok(x));

    while (s->length < limit && !fmpz_is_zero(x.b))
    {
        fmpz_poly_fit_length(s, s->length + 1);
        fmpz_fdiv_qr(s->coeffs + s->length, x.a, x.a, x.b);
        s->length++;
        fmpz_swap(x.a, x.b);
    }

    FLINT_ASSERT(!fmpz_is_zero(x.a));
    fmpz_swap(fmpq_numref(r), x.b);
    fmpz_swap(fmpq_denref(r), x.a);
    FLINT_ASSERT(fmpq_is_canonical(r));

    for (i = 0; i < s->length; i++)
        fmpz_swap(c + i, s->coeffs + i);

    return i;
}

/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_poly.h"

/*
    ball representing the closed interval

             a        a + da
    x = [ -------- , -------- ]
           b + db       b

    Most functions assume a>0, b>0, db>=0, da>=0 with the notable
    exception fmpq_ball_gt_one. The only arithmetic for which these balls are
    well-suited is the calculation of M(x) for M in GL2(ZZ)
*/
typedef struct {
    fmpz_t a, b, db, da;
} fmpq_ball_struct;

typedef fmpq_ball_struct fmpq_ball_t[1];


typedef struct {
    fmpz_t _11, _12, _21, _22;
    int det;    /* 0,1,or,-1: 0 -> dont know, 1 -> 1, -1 -> -1 */
} fmpz_mat22_struct;

typedef fmpz_mat22_struct fmpz_mat22_t[1];


void fmpz_mat22_print(const fmpz_mat22_t M)
{
    printf("mat[\n   ");
    fmpz_print(M->_11); printf(",\n   ");
    fmpz_print(M->_12); printf(",\n   ");
    fmpz_print(M->_21); printf(",\n   ");
    fmpz_print(M->_22); printf(",\n");
    printf("\n]\n");
}

void fmpz_mat22_init(fmpz_mat22_t M)
{
    fmpz_init(M->_11);
    fmpz_init(M->_12);
    fmpz_init(M->_21);
    fmpz_init(M->_22);
    M->det = 0;
}

void fmpz_mat22_clear(fmpz_mat22_t M)
{
    fmpz_clear(M->_11);
    fmpz_clear(M->_12);
    fmpz_clear(M->_21);
    fmpz_clear(M->_22);
}

void fmpz_mat22_one(fmpz_mat22_t M)
{
    fmpz_one(M->_11);
    fmpz_zero(M->_12);
    fmpz_zero(M->_21);
    fmpz_one(M->_22);
    M->det = 1;
}

int fmpz_mat22_is_one(fmpz_mat22_t M)
{
    return fmpz_is_one(M->_11)
        && fmpz_is_zero(M->_12)
        && fmpz_is_zero(M->_21)
        && fmpz_is_one(M->_22);
}

/* M = M.N */
void fmpz_mat22_rmul(fmpz_mat22_t M, const fmpz_mat22_t N)
{
    fmpz_t a, b, c, d;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);
    fmpz_mul(a, M->_11, N->_11); fmpz_addmul(a, M->_12, N->_21);
    fmpz_mul(b, M->_11, N->_12); fmpz_addmul(b, M->_12, N->_22);
    fmpz_mul(c, M->_21, N->_11); fmpz_addmul(c, M->_22, N->_21);
    fmpz_mul(d, M->_21, N->_12); fmpz_addmul(d, M->_22, N->_22);
    fmpz_swap(M->_11, a);
    fmpz_swap(M->_12, b);
    fmpz_swap(M->_21, c);
    fmpz_swap(M->_22, d);
    M->det *= N->det;
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(c);
    fmpz_clear(d);
}

/* M = M->[q 1; 1 0] */
void fmpz_mat22_rmul_elem(fmpz_mat22_t M, const fmpz_t q)
{
    fmpz_addmul(M->_12, M->_11, q);
    fmpz_addmul(M->_22, M->_21, q);
    fmpz_swap(M->_11, M->_12);
    fmpz_swap(M->_21, M->_22);
    M->det *= -1;
}


void fmpq_ball_init(fmpq_ball_t x)
{
    fmpz_init(x->a);
    fmpz_init(x->b);
    fmpz_init(x->db);
    fmpz_init(x->da);
}

void fmpq_ball_clear(fmpq_ball_t x)
{
    fmpz_clear(x->a);
    fmpz_clear(x->b);
    fmpz_clear(x->db);
    fmpz_clear(x->da);
}

void fmpq_ball_print(const fmpq_ball_t x)
{
    printf("Ball[");
    fmpz_print(x->a);
    flint_printf("/(");
    fmpz_print(x->b);
    flint_printf(" + ");
    fmpz_print(x->db);
    flint_printf("), (");
    fmpz_print(x->a);
    flint_printf(" + ");
    fmpz_print(x->da);
    flint_printf(")/");
    fmpz_print(x->b);
    flint_printf("]");
}

void fmpq_ball_swap(fmpq_ball_t x, fmpq_ball_t y)
{
   fmpq_ball_struct t = *x;
   *x = *y;
   *y = t;
}

/* is x canonical and bounded away from 1, i.e. 1 < x ? */
int fmpq_ball_gt_one(const fmpq_ball_t x)
{
    fmpz_t temp;
    int r;

    if (fmpz_sgn(x->a) <= 0)
        return 0;
    if (fmpz_sgn(x->b) <= 0)
        return 0;
    if (fmpz_sgn(x->db) < 0)
        return 0;
    if (fmpz_sgn(x->da) < 0)
        return 0;

    if (fmpz_cmp(x->b, x->a) >= 0)
        return 0;

    fmpz_init(temp);
    fmpz_add(temp, x->b, x->db);
    r = fmpz_cmp(temp, x->a) < 0;
    fmpz_clear(temp);
    return r;
}

/* does y contain x? */
int fmpq_ball_contains(const fmpq_ball_t y, const fmpq_ball_t x)
{
    int ans = 0;
    fmpz_t t1, t2, t3;

    fmpz_init(t1);
    fmpz_init(t2);
    fmpz_init(t3);

    fmpz_add(t3, x->b, x->db);
    fmpz_mul(t1, y->a, t3);
    fmpz_add(t3, y->b, y->db);
    fmpz_mul(t2, x->a, t3);
    if (fmpz_cmp(t1, t2) > 0)
        goto cleanup;

    fmpz_add(t3, x->a, x->da);
    fmpz_mul(t1, y->b, t3);
    fmpz_add(t3, y->a, y->da);
    fmpz_mul(t2, x->b, t3);
    if (fmpz_cmp(t1, t2) > 0)
        goto cleanup;

    ans = 1;

cleanup:

    fmpz_clear(t1);
    fmpz_clear(t2);
    fmpz_clear(t3);

    return ans;
}


/* y = m^-1(x) */
void fmpq_ball_apply_mat22_inv(
    fmpq_ball_t y,
    const fmpz_mat22_t M,
    const fmpq_ball_t x)
{
    if (M->det == 1)
    {
        fmpz_mul(y->a, x->a, M->_22);
        fmpz_submul(y->a, x->b, M->_12);
        fmpz_submul(y->a, x->db, M->_12);
        fmpz_mul(y->b, x->b, M->_11);
        fmpz_submul(y->b, x->a, M->_21);
        fmpz_submul(y->b, x->da, M->_21);
    }
    else
    {
        FLINT_ASSERT(M->det == -1);
        /* det = -1 swaps the endpoints */
        fmpz_mul(y->a, x->b, M->_12);
        fmpz_submul(y->a, x->a, M->_22);
        fmpz_submul(y->a, x->da, M->_22);
        fmpz_mul(y->b, x->a, M->_21);
        fmpz_submul(y->b, x->b, M->_11);
        fmpz_submul(y->b, x->db, M->_11);
    }
    fmpz_mul(y->db, x->db, M->_11);
    fmpz_addmul(y->db, x->da, M->_21);
    fmpz_mul(y->da, x->db, M->_12);
    fmpz_addmul(y->da, x->da, M->_22);
}

/* y = [q 1; 1 0]^-1(x) */
void fmpq_ball_apply_mat22_inv_elem(
    fmpq_ball_t y,
    const fmpz_t q,
    const fmpq_ball_t x)
{
    fmpz_set(y->a, x->b);
    fmpz_set(y->b, x->a);
    fmpz_submul(y->b, x->b, q);
    fmpz_submul(y->b, x->db, q);
    fmpz_mul(y->db, x->db, q);
    fmpz_add(y->db, y->db, x->da);
    fmpz_set(y->da, x->db);
}


/*
    The interface is entirely mutable and the terms are streamed to s
    The input x is assumed to be 1 < x and the output x' will be 1 < x'
                 1
    x = s1 + ----------
                         1
                 ... + -------
                              1
                        sn + ---
                              x'

    M will be multiplied on the right by the [si 1; 1 0] for the si appended.

    There is a ton of room for improvement in this function, but it already
    seems to be subquadratic, that is,

        time(cf(Fibonacci[2^n + 1]/Fibonacci[2^n])) << n^2

    Should try to approach the speed of gcd(x->a, x->b)

    Things to try:
    - the gauss iteration accumulates directly into M, which might not be that great
    - M is not necessarily needed when this function is called from the top level
    - after splitting, the fmpq_ball_apply_mat22_inv(y, N, x) can probably be
      optimized because we should know some of the leading bits of y
    - when k < 10000, say, it is probably better to use lemher, which is still
      quadratic, and let y be x >> (k - 2*FLINT_BITS) instead of x >> (k/2).
      The recursive call to fmpq_ball_get_cfrac(s, N, y, limit) will produce
      an N with entries probably bounded by FLINT_BITS.
*/
void fmpq_ball_get_cfrac(
    fmpz_poly_t s,
    fmpz_mat22_t M,
    fmpq_ball_t x,
    const slong limit)
{
    flint_bitcnt_t k;
    fmpz_t q, r;
    fmpq_ball_t y;
    fmpz_mat22_t N;

    FLINT_ASSERT(limit >= 0);

    fmpz_init(q);
    fmpz_init(r);
    fmpq_ball_init(y);
    fmpz_mat22_init(N);

again:

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_ball_gt_one(x));

    if (s->length >= limit)       
        goto cleanup;

    k = fmpz_bits(x->a);
    if (k >= FLINT_BITS - 1)
        goto split;

gauss:

    fmpz_fdiv_qr(q, r, x->a, x->b);
    FLINT_ASSERT(fmpz_sgn(q) > 0);

    fmpq_ball_apply_mat22_inv_elem(y, q, x);
    if (!fmpq_ball_gt_one(y))
        goto cleanup;

    fmpq_ball_swap(x, y);
    fmpz_mat22_rmul_elem(M, q);
    fmpz_poly_fit_length(s, s->length + 1);
    fmpz_swap(s->coeffs + s->length, q);

    s->length++;
    goto again;
    
split:

    fmpz_fdiv_q_2exp(y->a, x->a, k/2);
    fmpz_fdiv_q_2exp(y->b, x->b, k/2);
    fmpz_fdiv_q_2exp(y->db, x->db, k/2);
    fmpz_fdiv_q_2exp(y->da, x->da, k/2);
    fmpz_add_ui(y->db, y->db, 2); /* 1 for b and 1 for db */
    fmpz_add_ui(y->da, y->da, 2); /* ditto */
    if (!fmpq_ball_gt_one(y))
        goto gauss;
    FLINT_ASSERT(fmpq_ball_contains(y, x));
    fmpz_mat22_one(N);
    fmpq_ball_get_cfrac(s, N, y, limit);
    if (fmpz_mat22_is_one(N))
        goto gauss;
    fmpq_ball_apply_mat22_inv(y, N, x);
    fmpq_ball_swap(x, y);
    fmpz_mat22_rmul(M, N);
    goto again;

cleanup:

    fmpz_clear(q);
    fmpz_clear(r);
    fmpq_ball_clear(y);
    fmpz_mat22_clear(N);

    FLINT_ASSERT(fmpz_poly_length(s) <= limit);
    return;
}


slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t f, slong limit)
{
    int cmp;
    slong i;
    fmpq_ball_t x;
    fmpz_poly_t s;
    fmpz_mat22_t M;

    /* handle infinite input */
    if (fmpz_is_zero(fmpq_denref(f)))
    {
        fmpz_zero(fmpq_numref(rem));
        fmpz_one(fmpq_denref(rem));
        return 0;
    }

    FLINT_ASSERT(fmpq_is_canonical(f));

    if (limit <= 0)
    {
        if (fmpq_is_zero(f))
        {
            /* produce infinite output */
            fmpz_one(fmpq_numref(rem));
            fmpz_zero(fmpq_denref(rem));
        }
        else
        {
            fmpq_inv(rem, f);
        }
        return 0;
    }

    fmpz_mat22_init(M);
    fmpz_mat22_one(M);

    fmpq_ball_init(x);
    fmpz_set(x->a, fmpq_numref(f));
    fmpz_set(x->b, fmpq_denref(f));
    fmpz_zero(x->db);
    fmpz_zero(x->da);

    fmpz_poly_init(s);
    s->length = 0;

    cmp = fmpz_cmp(x->a, x->b);
    if (cmp > 0)
    {
        fmpq_ball_get_cfrac(s, M, x, limit);
    }
    else
    {
        fmpz_poly_fit_length(s, 1);
        if (cmp < 0 && fmpz_sgn(x->a) >= 0)
            fmpz_zero(s->coeffs + 0);
        else
            fmpz_fdiv_qr(s->coeffs + s->length, x->a, x->a, x->b);
        s->length = 1;
        fmpz_swap(x->a, x->b);

        if (!fmpz_is_zero(x->b))
            fmpq_ball_get_cfrac(s, M, x, limit);
    }

    FLINT_ASSERT(fmpz_is_zero(x->b) || fmpq_ball_gt_one(x));
    FLINT_ASSERT(fmpz_is_zero(x->db));
    FLINT_ASSERT(fmpz_is_zero(x->da));

    while (s->length < limit && !fmpz_is_zero(x->b))
    {
        fmpz_poly_fit_length(s, s->length + 1);
        fmpz_fdiv_qr(s->coeffs + s->length, x->a, x->a, x->b);
        s->length++;
        fmpz_swap(x->a, x->b);
    }

    FLINT_ASSERT(!fmpz_is_zero(x->a));
    fmpz_swap(fmpq_numref(rem), x->b);
    fmpz_swap(fmpq_denref(rem), x->a);
    FLINT_ASSERT(fmpq_is_canonical(rem));

    FLINT_ASSERT(s->length <= limit);
    for (i = 0; i < s->length; i++)
        fmpz_swap(c + i, s->coeffs + i);

    fmpz_mat22_clear(M);
    fmpq_ball_clear(x);
    fmpz_poly_clear(s);

    return i;
}



slong fmpq_get_cfrac_naive(fmpz * c, fmpq_t rem, const fmpq_t x, slong n)
{
    fmpz_t p, q;
    slong i;

    fmpz_init(p);
    fmpz_init(q);

    fmpz_set(p, fmpq_numref(x));
    fmpz_set(q, fmpq_denref(x));

    for (i = 0; i < n && !fmpz_is_zero(q); i++)
    {
        fmpz_fdiv_qr(c + i, p, p, q);
        fmpz_swap(p, q);
    }

    fmpz_set(fmpq_numref(rem), q);
    fmpz_set(fmpq_denref(rem), p);
    fmpq_canonicalise(rem);

    fmpz_clear(p);
    fmpz_clear(q);

    return i;
}

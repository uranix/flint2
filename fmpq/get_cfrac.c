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

#define SWORD_MAX ((WORD(1) << (FLINT_BITS - 2)) + ((WORD(1) << (FLINT_BITS - 2)) - 1))

/*
    Currently 2 types for balls because fmpq_gball_t is more general than
    fmpq_ball_t but slower
*/

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

/*
    ball representing the closed interval [left, right]
*/
typedef struct {
    fmpz_t left_num, left_den, right_num, right_den;
    int exact;
} fmpq_gball_struct;

typedef fmpq_gball_struct fmpq_gball_t[1];

typedef struct {
    fmpz_t _11, _12, _21, _22;
    int det;    /* 0,1,or,-1: 0 -> dont know, 1 -> 1, -1 -> -1 */
} fmpz_mat22_struct;

typedef fmpz_mat22_struct fmpz_mat22_t[1];

flint_bitcnt_t fmpz_mat22_bits(const fmpz_mat22_t N)
{
    flint_bitcnt_t b = fmpz_bits(N->_11);
    b = FLINT_MAX(b, fmpz_bits(N->_12));
    b = FLINT_MAX(b, fmpz_bits(N->_21));
    b = FLINT_MAX(b, fmpz_bits(N->_22));
    return b;
}

void fmpz_mat22_print(const fmpz_mat22_t M)
{
    printf("mat[\n   ");
    fmpz_print(M->_11); printf(",\n   ");
    fmpz_print(M->_12); printf(",\n   ");
    fmpz_print(M->_21); printf(",\n   ");
    fmpz_print(M->_22); printf("]\n");
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

/* M = M*[q 1; 1 0] */
void fmpz_mat22_rmul_elem(fmpz_mat22_t M, const fmpz_t q)
{
    fmpz_addmul(M->_12, M->_11, q);
    fmpz_addmul(M->_22, M->_21, q);
    fmpz_swap(M->_11, M->_12);
    fmpz_swap(M->_21, M->_22);
    M->det *= -1;
}

/* M = [q 1; 1 0]*M */
void fmpz_mat22_lmul_elem(fmpz_mat22_t M, const fmpz_t q)
{
    fmpz_addmul(M->_21, M->_11, q);
    fmpz_addmul(M->_22, M->_12, q);
    fmpz_swap(M->_11, M->_21);
    fmpz_swap(M->_12, M->_22);
    M->det *= -1;
}


void fmpq_ball_init(fmpq_ball_t x)
{
    fmpz_init(x->a);
    fmpz_init(x->b);
    fmpz_init(x->db);
    fmpz_init(x->da);
}

void fmpq_gball_init(fmpq_gball_t x)
{
    fmpz_init(x->left_num);
    fmpz_init(x->left_den);
    fmpz_init(x->right_num);
    fmpz_init(x->right_den);
    x->exact = 0;
}

void fmpq_ball_clear(fmpq_ball_t x)
{
    fmpz_clear(x->a);
    fmpz_clear(x->b);
    fmpz_clear(x->db);
    fmpz_clear(x->da);
}

void fmpq_gball_clear(fmpq_gball_t x)
{
    fmpz_clear(x->left_num);
    fmpz_clear(x->left_den);
    fmpz_clear(x->right_num);
    fmpz_clear(x->right_den);
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

void fmpq_gball_print(const fmpq_gball_t x)
{
    if (x->exact)
        printf("ExactBall[");
    else
        printf("Ball[");
    fmpz_print(x->left_num);
    flint_printf("/");
    fmpz_print(x->left_den);
    flint_printf(", ");
    fmpz_print(x->right_num);
    flint_printf("/");
    fmpz_print(x->right_den);
    flint_printf("]");
}

void fmpq_ball_swap(fmpq_ball_t x, fmpq_ball_t y)
{
   fmpq_ball_struct t = *x;
   *x = *y;
   *y = t;
}

void fmpq_gball_swap(fmpq_gball_t x, fmpq_gball_t y)
{
   fmpq_gball_struct t = *x;
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

int fmpq_gball_gt_one(const fmpq_gball_t x)
{
    if (fmpz_sgn(x->left_num) <= 0)
        return 0;
    if (fmpz_sgn(x->left_den) <= 0)
        return 0;
    if (fmpz_cmp(x->left_den, x->left_num) >= 0)
        return 0;

    if (x->exact)
        return 1;

    if (fmpz_sgn(x->right_num) <= 0)
        return 0;
    if (fmpz_sgn(x->right_den) <= 0)
        return 0;
    if (fmpz_cmp(x->right_den, x->right_num) >= 0)
        return 0;

    return 1;
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

int fmpq_gball_contains(const fmpq_gball_t y, const fmpq_gball_t x)
{
    if (_fmpq_cmp(y->left_num, y->left_den, x->left_num, x->left_den) > 0)
    {
        return 0;
    }

    if (_fmpq_cmp(x->exact ? x->left_num : x->right_num,
                  x->exact ? x->left_den : x->right_den,
                  y->exact ? y->left_num : y->right_num,
                  y->exact ? y->left_num : y->right_den) > 0)
    {
        return 0;
    }

    return 1;
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

void fmpq_gball_apply_mat22_inv(
    fmpq_gball_t y,
    const fmpz_mat22_t M,
    const fmpq_gball_t x)
{
    y->exact = x->exact;

    if (M->det == 1)
    {
        if (!x->exact)
        {
            fmpz_mul(y->right_num, x->right_num, M->_22);
            fmpz_submul(y->right_num, x->right_den, M->_12);
            fmpz_mul(y->right_den, x->right_den, M->_11);
            fmpz_submul(y->right_den, x->right_num, M->_21);
        }

        fmpz_mul(y->left_num, x->left_num, M->_22);
        fmpz_submul(y->left_num, x->left_den, M->_12);
        fmpz_mul(y->left_den, x->left_den, M->_11);
        fmpz_submul(y->left_den, x->left_num, M->_21);
    }
    else
    {
        FLINT_ASSERT(M->det == -1);
        /* det = -1 swaps the endpoints */

        if (x->exact)
        {
            fmpz_mul(y->left_num, x->left_den, M->_12);
            fmpz_submul(y->left_num, x->left_num, M->_22);
            fmpz_mul(y->left_den, x->left_num, M->_21);
            fmpz_submul(y->left_den, x->left_den, M->_11);
        }
        else
        {
            fmpz_mul(y->right_num, x->left_den, M->_12);
            fmpz_submul(y->right_num, x->left_num, M->_22);
            fmpz_mul(y->right_den, x->left_num, M->_21);
            fmpz_submul(y->right_den, x->left_den, M->_11);

            fmpz_mul(y->left_num, x->right_den, M->_12);
            fmpz_submul(y->left_num, x->right_num, M->_22);
            fmpz_mul(y->left_den, x->right_num, M->_21);
            fmpz_submul(y->left_den, x->right_den, M->_11);
        }
    }
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

void fmpq_gball_apply_mat22_inv_elem(
    fmpq_gball_t y,
    const fmpz_t q,
    const fmpq_gball_t x)
{
    y->exact = x->exact;

    if (x->exact)
    {
        fmpz_set(y->left_num, x->left_den);
        fmpz_set(y->left_den, x->left_num);
        fmpz_submul(y->left_den, x->left_den, q);
    }
    else
    {
        fmpz_set(y->right_num, x->left_den);
        fmpz_set(y->right_den, x->left_num);
        fmpz_submul(y->right_den, x->left_den, q);

        fmpz_set(y->left_num, x->right_den);
        fmpz_set(y->left_den, x->right_num);
        fmpz_submul(y->left_den, x->right_den, q);
    }
}


void fmpq_gball_apply_mat22_inv_elem2(
    fmpq_gball_t y,
    const fmpz_t q, fmpz_t r,
    const fmpq_gball_t x)
{
    y->exact = x->exact;

    if (x->exact)
    {
        fmpz_set(y->left_num, x->left_den);
/*
        fmpz_set(y->left_den, x->left_num);
        fmpz_submul(y->left_den, x->left_den, q);
        FLINT_ASSERT(fmpz_equal(y->left_den, r));
*/
        fmpz_swap(y->left_den, r);
    }
    else
    {
        fmpz_set(y->right_num, x->left_den);
/*
        fmpz_set(y->right_den, x->left_num);
        fmpz_submul(y->right_den, x->left_den, q);
        FLINT_ASSERT(fmpz_equal(y->right_den, r));
*/
        fmpz_swap(y->right_den, r);

        fmpz_set(y->left_num, x->right_den);
        fmpz_set(y->left_den, x->right_num);
        fmpz_submul(y->left_den, x->right_den, q);
    }
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
    - when k < 10000, say, it is probably better to use lehmer, which is still
      quadratic, and let y be x >> (k - 2*FLINT_BITS) instead of x >> (k/2).
      The recursive call to fmpq_ball_get_cfrac(s, N, y, limit) will produce
      an N with entries probably bounded by FLINT_BITS.
*/

void fmpq_ball_get_cfracOLD(
    fmpz_poly_t s,
    fmpz_mat22_t M, int needM,
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

    fmpz_mat22_one(M);

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
    if (needM)
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

    fmpq_ball_get_cfracOLD(s, N, 1, y, limit);
    if (fmpz_mat22_is_one(N))
        goto gauss;

    fmpq_ball_apply_mat22_inv(y, N, x);
    fmpq_ball_swap(x, y);
    if (needM)
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


void fmpq_ball_get_cfrac(
    fmpz_poly_t s,
    fmpz_mat22_t M, int needM,
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
    fmpz_mat22_one(M);

again:

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_ball_gt_one(x));

    if (s->length >= limit)
        goto cleanup;

    k = fmpz_bits(x->a);

    if (k > 4*FLINT_BITS)
    {
        if (k > 40*FLINT_BITS)
        {
            goto split;
        }
        else
        {
            goto lehmer;
        }
    }

gauss:

    if (s->length >= limit)
        goto cleanup;

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_ball_gt_one(x));

    fmpz_fdiv_qr(q, r, x->a, x->b);
    FLINT_ASSERT(fmpz_sgn(q) > 0);

    fmpq_ball_apply_mat22_inv_elem(y, q, x);
    if (!fmpq_ball_gt_one(y))
        goto cleanup;

    fmpq_ball_swap(x, y);
    if (needM)
        fmpz_mat22_rmul_elem(M, q);
    fmpz_poly_fit_length(s, s->length + 1);
    fmpz_swap(s->coeffs + s->length, q);
    s->length++;
    goto gauss;

lehmer:

    if (s->length >= limit)
        goto cleanup;

    k = fmpz_bits(x->a);
    if (k <= 4*FLINT_BITS)
        goto gauss;

    k -= 2*FLINT_BITS;

    fmpz_fdiv_q_2exp(y->a, x->a, k);
    fmpz_fdiv_q_2exp(y->b, x->b, k);
    fmpz_fdiv_q_2exp(y->db, x->db, k);
    fmpz_fdiv_q_2exp(y->da, x->da, k);
    fmpz_add_ui(y->db, y->db, 2); /* 1 for b and 1 for db */
    fmpz_add_ui(y->da, y->da, 2); /* ditto */
    if (!fmpq_ball_gt_one(y))
        goto single_step;

    fmpq_ball_get_cfrac(s, N, 1, y, limit);
    if (fmpz_mat22_is_one(N))
        goto single_step;

    fmpq_ball_apply_mat22_inv(y, N, x);
    if (needM)
        fmpz_mat22_rmul(M, N);
    fmpq_ball_swap(x, y);
    goto lehmer;

single_step:

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_ball_gt_one(x));

    fmpz_fdiv_qr(q, r, x->a, x->b);
    FLINT_ASSERT(fmpz_sgn(q) > 0);

    fmpq_ball_apply_mat22_inv_elem(y, q, x);
    if (!fmpq_ball_gt_one(y))
        goto cleanup;

    fmpq_ball_swap(x, y);
    if (needM)
        fmpz_mat22_rmul_elem(M, q);
    fmpz_poly_fit_length(s, s->length + 1);
    fmpz_swap(s->coeffs + s->length, q);
    s->length++;
    goto again;

split:

    k = fmpz_bits(x->a);
    k = k/2 + 2;
    k = FLINT_MAX(k, fmpz_bits(x->db));
    k -= 2;

    fmpz_fdiv_q_2exp(y->a, x->a, k);
    fmpz_fdiv_q_2exp(y->b, x->b, k);
    fmpz_fdiv_q_2exp(y->db, x->db, k);
    fmpz_fdiv_q_2exp(y->da, x->da, k);
    fmpz_add_ui(y->db, y->db, 2); /* 1 for b and 1 for db */
    fmpz_add_ui(y->da, y->da, 2); /* ditto */

    if (!fmpq_ball_gt_one(y))
        goto single_step;

    /*FLINT_ASSERT(fmpq_ball_contains(y, x));*/
    fmpq_ball_get_cfrac(s, N, 1, y, limit);
    if (fmpz_mat22_is_one(N))
        goto single_step;

    fmpq_ball_apply_mat22_inv(y, N, x);

    fmpz_mat22_rmul(M, N);

    fmpq_ball_swap(x, y);
    fmpq_ball_get_cfrac(s, N, needM, x, limit);

    if (needM)
        fmpz_mat22_rmul(M, N);

cleanup:

    fmpz_clear(q);
    fmpz_clear(r);
    fmpq_ball_clear(y);
    fmpz_mat22_clear(N);

    FLINT_ASSERT(fmpz_poly_length(s) <= limit);
    return;    
}


void fmpq_gball_chop(fmpq_gball_t y, const fmpq_gball_t x, ulong k)
{
    if (x->exact)
    {
        fmpz_fdiv_q_2exp(y->left_num, x->left_num, k);
        fmpz_fdiv_q_2exp(y->right_den, x->left_den, k);
        fmpz_add_ui(y->right_num, y->left_num, 1);
        fmpz_add_ui(y->left_den, y->right_den, 1);
    }
    else
    {
        fmpz_fdiv_q_2exp(y->left_num, x->left_num, k);
        fmpz_fdiv_q_2exp(y->left_den, x->left_den, k);
        fmpz_add_ui(y->left_den, y->left_den, 1);
        fmpz_fdiv_q_2exp(y->right_num, x->right_num, k);
        fmpz_fdiv_q_2exp(y->right_den, x->right_den, k);
        fmpz_add_ui(y->right_num, y->right_num, 1);
    }
    y->exact = 0;
}


void fmpq_gball_get_cfrac(
    fmpz_poly_t s,
    fmpz_mat22_t M, int needM,
    fmpq_gball_t x,
    const slong limit)
{
    flint_bitcnt_t k;
    fmpz_t q, r;
    fmpq_gball_t y;
    fmpz_mat22_t N;

    FLINT_ASSERT(limit >= 0);

    fmpz_init(q);
    fmpz_init(r);
    fmpq_gball_init(y);
    fmpz_mat22_init(N);

    fmpz_mat22_one(M);

    FLINT_ASSERT(fmpq_gball_gt_one(x));

again:

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(x->exact || _fmpq_cmp(x->left_num, x->left_den,
                                       x->right_num, x->right_den) <= 0);
    FLINT_ASSERT(fmpq_gball_gt_one(x));

    if (s->length >= limit)
        goto cleanup;

    k = fmpz_bits(x->left_num);

    if (k > 4*FLINT_BITS)
    {
        if (k > 100*FLINT_BITS)
        {
            goto split;
        }
        else
        {
/*
            flint_printf("lehmer bits: %wu %wu %wu %wu\n", fmpz_bits(x->left_num),
                                                           fmpz_bits(x->left_den),
                                                           fmpz_bits(x->right_num),
                                                           fmpz_bits(x->right_den));
*/
            goto lehmer;
        }
    }

/*
            flint_printf("gauss bits: %wu %wu %wu %wu\n", fmpz_bits(x->left_num),
                                                           fmpz_bits(x->left_den),
                                                           fmpz_bits(x->right_num),
                                                           fmpz_bits(x->right_den));
*/


gauss:

    if (s->length >= limit)
        goto cleanup;

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_gball_gt_one(x));

    fmpz_fdiv_qr(q, r, x->left_num, x->left_den);
    FLINT_ASSERT(fmpz_sgn(q) > 0);

    fmpq_gball_apply_mat22_inv_elem2(y, q, r, x);
    if (!fmpq_gball_gt_one(y))
        goto cleanup;

    fmpq_gball_swap(x, y);

    if (needM)
        fmpz_mat22_rmul_elem(M, q);

    fmpz_poly_fit_length(s, s->length + 1);
    fmpz_swap(s->coeffs + s->length, q);
    s->length++;
    goto gauss;

lehmer:

    FLINT_ASSERT(fmpq_gball_gt_one(x));

    if (s->length >= limit)
        goto cleanup;

    k = fmpz_bits(x->left_num);
    if (k <= 4*FLINT_BITS)
        goto gauss;

    k -= 2*FLINT_BITS;
    fmpq_gball_chop(y, x, k);

    if (!fmpq_gball_gt_one(y))
        goto single_step;

    fmpq_gball_get_cfrac(s, N, 1, y, limit);
    if (fmpz_mat22_is_one(N))
        goto single_step;

    fmpq_gball_apply_mat22_inv(y, N, x);

    if (needM)
        fmpz_mat22_rmul(M, N);

    fmpq_gball_swap(x, y);
    FLINT_ASSERT(fmpq_gball_gt_one(x));

    goto lehmer;

single_step:

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_gball_gt_one(x));

    fmpz_fdiv_qr(q, r, x->left_num, x->left_den);
    FLINT_ASSERT(fmpz_sgn(q) > 0);

    fmpq_gball_apply_mat22_inv_elem2(y, q, r, x);
    if (!fmpq_gball_gt_one(y))
        goto cleanup;

    fmpq_gball_swap(x, y);

    if (needM)
        fmpz_mat22_rmul_elem(M, q);

    fmpz_poly_fit_length(s, s->length + 1);
    fmpz_swap(s->coeffs + s->length, q);
    s->length++;
    goto again;

split:

    FLINT_ASSERT(fmpq_gball_gt_one(x));

    k = fmpz_bits(x->left_num);
    k = k/2 + 2;

    fmpq_gball_chop(y, x, k);
    if (!fmpq_gball_gt_one(y))
        goto single_step;

/*    FLINT_ASSERT(fmpq_gball_contains(y, x));*/
    fmpq_gball_get_cfrac(s, N, 1, y, limit);
    if (fmpz_mat22_is_one(N))
        goto single_step;

    fmpq_gball_apply_mat22_inv(y, N, x);

    if (needM)
        fmpz_mat22_rmul(M, N);

    fmpq_gball_swap(x, y);
    FLINT_ASSERT(fmpq_gball_gt_one(x));

    fmpq_gball_get_cfrac(s, N, needM, x, limit);

    if (needM)
        fmpz_mat22_rmul(M, N);

cleanup:

    fmpz_clear(q);
    fmpz_clear(r);
    fmpq_gball_clear(y);
    fmpz_mat22_clear(N);

    FLINT_ASSERT(fmpz_poly_length(s) <= limit);
    return;
}




void fmpq_simplest_between(fmpq_t mid, const fmpq_t l, const fmpq_t r)
{
    if (fmpq_cmp(l, r) <= 0)
    {
        _fmpq_simplest_between(fmpq_numref(mid), fmpq_denref(mid),
                               fmpq_numref(l), fmpq_denref(l),
                               fmpq_numref(r), fmpq_denref(r));
    }
    else
    {
        _fmpq_simplest_between(fmpq_numref(mid), fmpq_denref(mid),
                               fmpq_numref(r), fmpq_denref(r),
                               fmpq_numref(l), fmpq_denref(l));
    }
}

void _fmpq_simplest_between(
    fmpz_t mid_num, fmpz_t mid_den,
    const fmpz_t l_num, const fmpz_t l_den,
    const fmpz_t r_num, const fmpz_t r_den)
{
    fmpz_t q;
    fmpz_poly_t s;
    fmpz_mat22_t M;
    fmpq_gball_t x;

    FLINT_ASSERT(fmpz_sgn(l_den) > 0);
    FLINT_ASSERT(fmpz_sgn(r_den) > 0);
    FLINT_ASSERT(_fmpq_cmp(l_num, l_den, r_num, r_den) <= 0);

    fmpz_init(q);
    fmpz_poly_init(s);
    fmpq_gball_init(x);
    fmpz_mat22_init(M);

    fmpz_mat22_one(M);

    fmpz_set(x->left_num, l_num);
    fmpz_set(x->left_den, l_den);
    fmpz_set(x->right_num, r_num);
    fmpz_set(x->right_den, r_den);

    if (fmpz_cmp(x->left_num, x->left_den) > 0)
    {
        /* 1 < x */
        fmpq_gball_get_cfrac(s, M, 1, x, SWORD_MAX);
    }
    else if (fmpz_sgn(x->left_num) > 0
                 && fmpz_cmp(x->right_num, x->right_den) < 0)
    {
        /* 0 < x < 1 */
        fmpz_swap(x->left_den, x->right_num);
        fmpz_swap(x->left_num, x->right_den);
        fmpq_gball_get_cfrac(s, M, 1, x, SWORD_MAX);
        fmpz_zero(q);
        fmpz_mat22_lmul_elem(M, q);
    }
    else
    {
        fmpq_gball_t y;
        fmpq_gball_init(y);
        fmpz_fdiv_q(q, x->left_num, x->left_den);
        fmpq_gball_apply_mat22_inv_elem(y, q, x);
        if (fmpq_gball_gt_one(y))
        {
            fmpq_gball_swap(x, y);
            fmpq_gball_get_cfrac(s, M, 1, x, SWORD_MAX);
            fmpz_mat22_lmul_elem(M, q);
        }
        fmpq_gball_clear(y);
    }

    fmpz_cdiv_q(q, x->left_num, x->left_den);

#if WANT_ASSERT
    {
        fmpz_t one;
        fmpz_init_set_ui(one, 1);
        FLINT_ASSERT(_fmpq_cmp(x->right_num, x->right_den, q, one) >= 0);
        fmpz_clear(one);
    }
#endif

    FLINT_ASSERT(M->det == 1 || M->det == -1);
    fmpz_swap(mid_num, M->_12);
    fmpz_addmul(mid_num, M->_11, q);
    fmpz_swap(mid_den, M->_22);
    fmpz_addmul(mid_den, M->_21, q);

    fmpz_clear(q);
    fmpz_poly_clear(s);
    fmpq_gball_clear(x);
    fmpz_mat22_clear(M);

    FLINT_ASSERT(_fmpq_is_canonical(mid_num, mid_den));
    return;
}




slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t f, slong limit)
{
    slong i;
    int cmp, num_sign, den_sign;
    fmpq_gball_t x;
    fmpz_poly_t s;
    fmpz_mat22_t M;
#if WANT_ASSERT
    int input_is_canonical;
#endif

    num_sign = fmpz_sgn(fmpq_denref(f));
    den_sign = fmpz_sgn(fmpq_denref(f));

    if (limit <= 0 || den_sign == 0)
    {
        if (num_sign < 0)
        {
            fmpz_neg(fmpq_numref(rem), fmpq_numref(f));
            fmpz_neg(fmpq_denref(rem), fmpq_denref(f));
        }
        else
        {
            fmpz_set(fmpq_numref(rem), fmpq_numref(f));
            fmpz_set(fmpq_denref(rem), fmpq_denref(f));
        }
        fmpz_swap(fmpq_numref(rem), fmpq_denref(rem));
        return 0;
    }

#if WANT_ASSERT
    input_is_canonical = fmpq_is_canonical(f);
#endif

    fmpz_mat22_init(M);
    fmpz_mat22_one(M);

    fmpq_gball_init(x);
    if (den_sign > 0)
    {
        fmpz_set(x->left_num, fmpq_numref(f));
        fmpz_set(x->left_den, fmpq_denref(f));
    }
    else
    {
        fmpz_neg(x->left_num, fmpq_numref(f));
        fmpz_neg(x->left_den, fmpq_denref(f));
    }
    x->exact = 1;

    fmpz_poly_init(s);
    s->length = 0;

    cmp = fmpz_cmp(x->left_num, x->left_den);
    if (cmp > 0)
    {
        fmpq_gball_get_cfrac(s, M, 0, x, limit);
    }
    else
    {
        fmpz_poly_fit_length(s, 1);
        if (cmp < 0 && fmpz_sgn(x->left_num) >= 0)
        {
            fmpz_zero(s->coeffs + 0);
        }
        else
        {
            fmpz_fdiv_qr(s->coeffs + s->length, x->left_num, x->left_num, x->left_den);
        }
        s->length = 1;
        fmpz_swap(x->left_num, x->left_den);

        if (!fmpz_is_zero(x->left_den))
            fmpq_gball_get_cfrac(s, M, 0, x, limit);
    }

    FLINT_ASSERT(fmpz_is_zero(x->left_den) || fmpq_gball_gt_one(x));
    FLINT_ASSERT(x->exact);

    while (s->length < limit && !fmpz_is_zero(x->left_den))
    {
        fmpz_poly_fit_length(s, s->length + 1);
        fmpz_fdiv_qr(s->coeffs + s->length, x->left_num, x->left_num, x->left_den);
        s->length++;
        fmpz_swap(x->left_num, x->left_den);
    }

    /* write remainder */
    FLINT_ASSERT(!fmpz_is_zero(x->left_num));
    fmpz_swap(fmpq_numref(rem), x->left_den);
    fmpz_swap(fmpq_denref(rem), x->left_num);
    FLINT_ASSERT(!input_is_canonical || fmpq_is_canonical(rem));

    /* write terms */
    FLINT_ASSERT(s->length <= limit);
    for (i = 0; i < s->length; i++)
        fmpz_swap(c + i, s->coeffs + i);

    fmpz_mat22_clear(M);
    fmpq_gball_clear(x);
    fmpz_poly_clear(s);

    return i;
}

#if 0
slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t f, slong limit)
{
    slong i;
    int cmp, num_sign, den_sign;
    fmpq_ball_t x;
    fmpz_poly_t s;
    fmpz_mat22_t M;
#if WANT_ASSERT
    int input_is_canonical;
#endif

    num_sign = fmpz_sgn(fmpq_denref(f));
    den_sign = fmpz_sgn(fmpq_denref(f));

    if (limit <= 0 || den_sign == 0)
    {
        if (num_sign < 0)
        {
            fmpz_neg(fmpq_numref(rem), fmpq_numref(f));
            fmpz_neg(fmpq_denref(rem), fmpq_denref(f));
        }
        else
        {
            fmpz_set(fmpq_numref(rem), fmpq_numref(f));
            fmpz_set(fmpq_denref(rem), fmpq_denref(f));
        }
        fmpz_swap(fmpq_numref(rem), fmpq_denref(rem));
        return 0;
    }

#if WANT_ASSERT
    input_is_canonical = fmpq_is_canonical(f);
#endif

    fmpz_mat22_init(M);
    fmpz_mat22_one(M);

    fmpq_ball_init(x);
    if (den_sign > 0)
    {
        fmpz_set(x->a, fmpq_numref(f));
        fmpz_set(x->b, fmpq_denref(f));
    }
    else
    {
        fmpz_neg(x->a, fmpq_numref(f));
        fmpz_neg(x->b, fmpq_denref(f));
    }
    fmpz_zero(x->db);
    fmpz_zero(x->da);

    fmpz_poly_init(s);
    s->length = 0;

    cmp = fmpz_cmp(x->a, x->b);
    if (cmp > 0)
    {
        fmpq_ball_get_cfrac(s, M, 0, x, limit);
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
            fmpq_ball_get_cfrac(s, M, 0, x, limit);
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

    /* write remainder */
    FLINT_ASSERT(!fmpz_is_zero(x->a));
    fmpz_swap(fmpq_numref(rem), x->b);
    fmpz_swap(fmpq_denref(rem), x->a);
    FLINT_ASSERT(!input_is_canonical || fmpq_is_canonical(rem));

    /* write terms */
    FLINT_ASSERT(s->length <= limit);
    for (i = 0; i < s->length; i++)
        fmpz_swap(c + i, s->coeffs + i);

    fmpz_mat22_clear(M);
    fmpq_ball_clear(x);
    fmpz_poly_clear(s);

    return i;
}
#endif


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

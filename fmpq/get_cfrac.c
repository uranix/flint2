/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "longlong.h"
#include "fmpq.h"
#include "fmpz_poly.h"

#define MPN_EXTRACT_NUMB(count, xh, xl)				\
  ((((xh) << ((count) - GMP_NAIL_BITS)) & GMP_NUMB_MASK) |	\
   ((xl) >> (GMP_LIMB_BITS - (count))))

/* Realloc for an mpz_t WHAT if it has less than NEEDED limbs.  */
#define MPZ_REALLOC(z,n) ((n) > ALLOC(z)			\
			  ? (mp_ptr) _mpz_realloc(z,n)			\
			  : PTR(z))
#define MPZ_NEWALLOC(z,n) ((n) > ALLOC(z)			\
			   ? (mp_ptr) _mpz_newalloc(z,n)		\
			   : PTR(z))

/* Field access macros.  */
#define SIZ(x) ((x)->_mp_size)
#define ABSIZ(x) ABS (SIZ (x))
#define PTR(x) ((x)->_mp_d)
#define EXP(x) ((x)->_mp_exp)
#define PREC(x) ((x)->_mp_prec)
#define ALLOC(x) ((x)->_mp_alloc)

#define MP_LIMB_T_SWAP(x, y)						\
  do {									\
    mp_limb_t __mp_limb_t_swap__tmp = (x);				\
    (x) = (y);								\
    (y) = __mp_limb_t_swap__tmp;					\
  } while (0)
#define MP_SIZE_T_SWAP(x, y)						\
  do {									\
    mp_size_t __mp_size_t_swap__tmp = (x);				\
    (x) = (y);								\
    (y) = __mp_size_t_swap__tmp;					\
  } while (0)


#define MP_SRCPTR_SWAP(x, y)						\
  do {									\
    mp_srcptr __mp_srcptr_swap__tmp = (x);				\
    (x) = (y);								\
    (y) = __mp_srcptr_swap__tmp;					\
  } while (0)

#define MPN_PTR_SWAP(xp,xs, yp,ys)					\
  do {									\
    MP_PTR_SWAP (xp, yp);						\
    MP_SIZE_T_SWAP (xs, ys);						\
  } while(0)
#define MPN_SRCPTR_SWAP(xp,xs, yp,ys)					\
  do {									\
    MP_SRCPTR_SWAP (xp, yp);						\
    MP_SIZE_T_SWAP (xs, ys);						\
  } while(0)

#define MPZ_PTR_SWAP(x, y)						\
  do {									\
    mpz_ptr __mpz_ptr_swap__tmp = (x);					\
    (x) = (y);								\
    (y) = __mpz_ptr_swap__tmp;						\
  } while (0)
#define MPZ_SRCPTR_SWAP(x, y)						\
  do {									\
    mpz_srcptr __mpz_srcptr_swap__tmp = (x);				\
    (x) = (y);								\
    (y) = __mpz_srcptr_swap__tmp;					\
  } while (0)


typedef struct {
    fmpz_t left_num, left_den, right_num, right_den;
    int exact;
} fmpq_ball_struct;

typedef fmpq_ball_struct fmpq_ball_t[1];

typedef struct {
    fmpz_t _11, _12, _21, _22;
    int det;    /* 0,1,or,-1: 0 -> dont know, 1 -> 1, -1 -> -1 */
} fmpz_mat22_struct;

typedef fmpz_mat22_struct fmpz_mat22_t[1];

typedef struct {
    mp_limb_t _11, _12, _21, _22;
    int det;
} ui_mat22_struct;

typedef ui_mat22_struct ui_mat22_t[1];

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

void ui_mat22_print(const ui_mat22_t M)
{
    flint_printf("mat[{{%wu, %wu}, {%wu, %wu}}]", M->_11, M->_12, M->_21, M->_22);
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

void fmpz_mat22_rmul_ui(fmpz_mat22_t M, const ui_mat22_t N)
{
    fmpz_t a, b, c, d;
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);
    fmpz_init(d);
    fmpz_mul_ui(a, M->_11, N->_11); fmpz_addmul_ui(a, M->_12, N->_21);
    fmpz_mul_ui(b, M->_11, N->_12); fmpz_addmul_ui(b, M->_12, N->_22);
    fmpz_mul_ui(c, M->_21, N->_11); fmpz_addmul_ui(c, M->_22, N->_21);
    fmpz_mul_ui(d, M->_21, N->_12); fmpz_addmul_ui(d, M->_22, N->_22);
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
    fmpz_init(x->left_num);
    fmpz_init(x->left_den);
    fmpz_init(x->right_num);
    fmpz_init(x->right_den);
    x->exact = 0;
}

void fmpq_ball_clear(fmpq_ball_t x)
{
    fmpz_clear(x->left_num);
    fmpz_clear(x->left_den);
    fmpz_clear(x->right_num);
    fmpz_clear(x->right_den);
}

void fmpq_ball_print(const fmpq_ball_t x)
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

/* is x canonical and bounded away from 1, i.e. 1 < x ? */
int fmpq_ball_gt_one(const fmpq_ball_t x)
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


void fmpq_ball_apply_mat22_inv_elem2(
    fmpq_ball_t y,
    const fmpz_t q, fmpz_t r,
    const fmpq_ball_t x)
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


void fmpq_ball_chop(fmpq_ball_t y, const fmpq_ball_t x, ulong k)
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


/*
    assume  2^(FLINT_BITS) <= a < 2^(2*FLINT_BITS)
            2^(FLINT_BITS) <= b < 2^(2*FLINT_BITS)
            0 < b <= a
    produce quotients s1, ..., sn, and M = prod_i [qi 1; 1 0]  with
        bits(M) <= FLINT_BITS such that s1, ..., sn is a prefix of the cfrac of
        any element of the open interval (a/(b+1), (a+1)/b)
    return n
*/
static slong _red2(
    mp_limb_t * s,
    mp_limb_t A1, mp_limb_t A0,
    mp_limb_t B1, mp_limb_t B0,
    ui_mat22_t M)
{
    int i;
    slong written = 0;
    mp_limb_t d0, d1;
    mp_limb_t t0, t1, t2, r0, r1;
    int det = 1;
    mp_limb_t m11 = 1;
    mp_limb_t m12 = 0;
    mp_limb_t m21 = 0;
    mp_limb_t m22 = 1;
    mp_limb_t a1 = A1;
    mp_limb_t a0 = A0;
    mp_limb_t b1 = B1;
    mp_limb_t b0 = B0;
    mp_limb_t q;

    FLINT_ASSERT(a1 != 0);
    FLINT_ASSERT(b1 < a1 || (b1 == a1 && b0 <= a0));

    if (b1 == 0 || !(b1 < a1 || (b1 == a1 && b0 < a0)))
        goto done;

    while (1)
    {
        FLINT_ASSERT(a1 != 0);
        FLINT_ASSERT(b1 != 0);
        FLINT_ASSERT(b1 <= a1);
        FLINT_ASSERT(b1 < a1 || (b1 == a1 && b0 < a0));
        q = 1;

        sub_ddmmss(r1,r0, a1,a0, b1,b0);

        for (i = 2; i <= 4; i++)
        {
            sub_dddmmmsss(t2,t1,t0, 0,r1,r0, 0,b1,b0);
            if (t2 != 0)
                goto quotient_found;
            q += 1;
            r0 = t0;
            r1 = t1;
        }

        if (r1 != 0)
        {
            int ncnt, dcnt;
            mp_limb_t Q = 0;

            count_leading_zeros(ncnt, r1);
            count_leading_zeros(dcnt, b1);
            dcnt -= ncnt;

            d1 = (b1 << dcnt) | ((-(ulong)(dcnt > 0)) & (b0 >> (FLINT_BITS - dcnt)));
            d0 = b0 << dcnt;

            while (dcnt >= 0)
            {
                sub_dddmmmsss(t2,t1,t0, 0,r1,r0, 0,d1,d0);
                Q = 2*Q + 1 + t2;
                r1 = t2 ? r1 : t1;
                r0 = t2 ? r0 : t0;
                d0 = (d1 << (FLINT_BITS - 1)) | (d0 >> 1);
                d1 = d1 >> 1;
                dcnt--;
            }

            q += Q;
        }

quotient_found:

        FLINT_ASSERT(r1 < b1 || (r1 == b1 && r0 < b0));

        t1 = m12 + q*m11;
        t2 = m22 + q*m21;

        if (r1 == 0)
            break;

        a0 = b0;
        a1 = b1;
        b0 = r0;
        b1 = r1;

        m12 = m11;
        m22 = m21;
        m11 = t1;
        m21 = t2;
        det *= -1;

        s[written] = q;
        written++;
        FLINT_ASSERT(written <= 2*FLINT_BITS);
    }

    FLINT_ASSERT(a1 != 0);
    FLINT_ASSERT(b1 <= a1);
    FLINT_ASSERT(b1 < a1 || (b1 == a1 && b0 < a0));

    /*
        if det 1
            big {A + 1, B} = {{m11, m12}, {m21, m22}} . {a + m22, b - m21} big
             sm {A, B + 1} = {{m11, m12}, {m21, m22}} . {a - m12, b + m11} sm
            need  b >= m21
                  a - m12 >= b + m11  i.e.  a - b >= m11 + m12

        else det -1
            big {A + 1, B} = {{m11, m12}, {m21, m22}} . {a - m22, b + m21} sm
             sm {A, B + 1} = {{m11, m12}, {m21, m22}} . {a + m12, b - m11} big
            need  b >= m11
                  a - m22 >= b + m21  i.e.  a - b >= m21 + m22
    */

    sub_ddmmss(d1,d0, a1,a0, b1,b0);
    if (det == 1)
    {
        if (b1 == 0 && b0 < m21)
            goto fix;
        add_ssaaaa(t1,t0, 0,m11, 0,m12);
    }
    else
    {
        if (b1 == 0 && b0 < m11)
            goto fix;        
        add_ssaaaa(t1,t0, 0,m21, 0,m22);
    }
    if (d1 < t1 || (d1 == t1 && d0 < t0))
        goto fix;

fixed:

#if WANT_ASSERT

    /* should be fixed */

    sub_ddmmss(d1,d0, a1,a0, b1,b0);
    if (det == 1)
    {
        if (b1 == 0 && b0 < m21)
            FLINT_ASSERT(0);
        add_ssaaaa(t1,t0, 0,m11, 0,m12);
    }
    else
    {
        if (b1 == 0 && b0 < m11)
            FLINT_ASSERT(0);
        add_ssaaaa(t1,t0, 0,m21, 0,m22);
    }
    if (d1 < t1 || (d1 == t1 && d0 < t0))
        FLINT_ASSERT(0);

    /* should have {A, B} == {{m11, m12}, {m21, m22}} . {a, b} */

    umul_ppmm(t1,t0, a0, m11); t1 += a1*m11;
    umul_ppmm(d1,d0, b0, m12); d1 += b1*m12;
    add_sssaaaaaa(t2,t1,t0, 0,t1,t0, 0,d1,d0);
    FLINT_ASSERT(t0 == A0 && t1 == A1 && t2 == 0);

    umul_ppmm(t1,t0, a0, m21); t1 += a1*m21;
    umul_ppmm(d1,d0, b0, m22); d1 += b1*m22;
    add_sssaaaaaa(t2,t1,t0, 0,t1,t0, 0,d1,d0);
    FLINT_ASSERT(t0 == B0 && t1 == B1 && t2 == 0);

    /* should have det = Det[{{m11, m12}, {m21, m22}}] */

    umul_ppmm(t1,t0, m11, m22);
    umul_ppmm(d1,d0, m12, m21);
    sub_ddmmss(t1,t0, t1,t0, d1,d0);
    FLINT_ASSERT(t1 == FLINT_SIGN_EXT(t0));
    FLINT_ASSERT(t0 == det);

#endif

done:

    M->_11 = m11;
    M->_12 = m12;
    M->_21 = m21;
    M->_22 = m22;
    M->det = det;

    return written;

fix:

    written--;
    FLINT_ASSERT(written >= 0);

    q = s[written];

    t1 = m11 - q*m12;
    t2 = m21 - q*m22;
    m11 = m12;
    m21 = m22;
    m12 = t1;
    m22 = t2;
    det *= -1;

#if WANT_ASSERT
    umul_ppmm(t1,t0, a0, q); t1 += a1*q;
    add_ssaaaa(t1,t0, t1,t0, b1,b0);
    b0 = a0;
    b1 = a1;
    a0 = t0;
    a1 = t1;
#endif

    goto fixed;

}

static void my_mpz_swap(mpz_ptr u, mpz_ptr v)
{
  MP_SIZE_T_SWAP(ALLOC(u), ALLOC(v));
  MP_SIZE_T_SWAP(SIZ(u), SIZ(v));
  MP_PTR_SWAP(PTR(v), PTR(u));
}

/* y = a1*x1 - a2*x2 */
static mp_size_t _msub(mp_ptr y, mp_limb_t a1, mp_ptr x1
                               , mp_limb_t a2, mp_ptr x2, mp_size_t n)
{
    mp_limb_t h0, h1;

    h0 =    mpn_mul_1(y, x1, n, a1);
    h1 = mpn_submul_1(y, x2, n, a2);

    if (h1 != h0)
        return -1;

    while (n > 0 && y[n - 1] == 0)
        n--;

    return n;
}

static void fmpq_ball_get_cfrac_lehmer_exact(
    fmpz_poly_t s, const slong limit,
    fmpz_mat22_t M, int needM,
    mpz_ptr xn,
    mpz_ptr xd,
    mpz_ptr yn,
    mpz_ptr yd)
{
    mp_limb_t s_temp[2*FLINT_BITS];
    slong i, written;
    mp_size_t xn_len;
    mp_size_t xd_len;
    mp_size_t yn_len;
    mp_size_t yd_len;
    mp_ptr xn_ptr;
    mp_ptr xd_ptr;
    mp_ptr yn_ptr;
    mp_ptr yd_ptr;
    ui_mat22_t m;
    mp_limb_t A0, A1, B0, B1;
    mp_size_t n;

    fmpz_mat22_one(M);

    /* fit everything to xn_len */
    n = SIZ(xn);

    MPZ_REALLOC(xd, n);
    MPZ_REALLOC(yn, n);
    MPZ_REALLOC(yd, n);

again:

    xn_len = SIZ(xn);
    xd_len = SIZ(xd);

    xn_ptr = PTR(xn);
    xd_ptr = PTR(xd);
    yn_ptr = PTR(yn);
    yd_ptr = PTR(yd);

    /* supposed xn > xd > 0 */
    FLINT_ASSERT(xn_len >= xd_len);
    FLINT_ASSERT(xd_len > 0);

    FLINT_ASSERT(xn_ptr[xn_len-1] != 0);
    FLINT_ASSERT(xd_ptr[xd_len-1] != 0);

    n = xn_len;
    if (n < 3)
        goto cleanup;

    if (n != xd_len && n != xd_len + 1)
        goto cleanup;

    if (n == xd_len + 1)
        xd_ptr[n-1] = 0;

    if ((slong)(xn_ptr[n-1]) < 0)
    {
        A1 = xn_ptr[n-1];
        A0 = xn_ptr[n-2];
        B1 = xd_ptr[n-1];
        B0 = xd_ptr[n-2];
    }
    else
    {
        int shift;
        count_leading_zeros(shift, xn_ptr[n - 1]);
        A1 = MPN_EXTRACT_NUMB(shift, xn_ptr[n-1], xn_ptr[n-2]);
        A0 = MPN_EXTRACT_NUMB(shift, xn_ptr[n-2], xn_ptr[n-3]);
        B1 = MPN_EXTRACT_NUMB(shift, xd_ptr[n-1], xd_ptr[n-2]);
        B0 = MPN_EXTRACT_NUMB(shift, xd_ptr[n-2], xd_ptr[n-3]);
    }

    written = _red2(s_temp, A1, A0, B1, B0, m);
    if (written <= 0 || s->length + written > limit)
        goto cleanup;

    if (m->det == 1)
    {
        yn_len = _msub(yn_ptr, m->_22, xn_ptr, m->_12, xd_ptr, n);
        if (yn_len <= 0)
            goto cleanup;

        yd_len = _msub(yd_ptr, m->_11, xd_ptr, m->_21, xn_ptr, n);
        if (yd_len <= 0)
            goto cleanup;
    }
    else
    {
        yn_len = _msub(yn_ptr, m->_12, xd_ptr, m->_22, xn_ptr, n);
        if (yn_len <= 0)
            goto cleanup;

        yd_len = _msub(yd_ptr, m->_21, xn_ptr, m->_11, xd_ptr, n);
        if (yd_len <= 0)
            goto cleanup;
    }

    SIZ(yn) = yn_len;
    SIZ(yd) = yd_len;

    if (needM)
        fmpz_mat22_rmul_ui(M, m);

    fmpz_poly_fit_length(s, s->length + written);
    for (i = 0; i < written; i++)
        fmpz_set_ui(s->coeffs + s->length + i, s_temp[i]);
    s->length += written;
    FLINT_ASSERT(s->length <= limit);

    my_mpz_swap(xn, yn);
    my_mpz_swap(xd, yd);
    goto again;

cleanup:

    /* never wrong */
    SIZ(yn) = 0;
    SIZ(yd) = 0;

    return;

}


static void fmpq_ball_get_cfrac_lehmer_inexact(
    fmpz_poly_t s, const slong limit,
    fmpz_mat22_t M, int needM,
    mpz_ptr xln,
    mpz_ptr xld,
    mpz_ptr xrn,
    mpz_ptr xrd,
    mpz_ptr yln,
    mpz_ptr yld,
    mpz_ptr yrn,
    mpz_ptr yrd)
{
    mp_limb_t s_temp[2*FLINT_BITS];
    slong i, written;
    mp_size_t xln_len;
    mp_size_t xld_len;
    mp_size_t xrn_len;
    mp_size_t xrd_len;
    mp_size_t yln_len;
    mp_size_t yld_len;
    mp_size_t yrn_len;
    mp_size_t yrd_len;
    mp_ptr xln_ptr;
    mp_ptr xld_ptr;
    mp_ptr xrn_ptr;
    mp_ptr xrd_ptr;
    mp_ptr yln_ptr;
    mp_ptr yld_ptr;
    mp_ptr yrn_ptr;
    mp_ptr yrd_ptr;
    ui_mat22_t m;
    mp_limb_t A0, A1, B0, B1;
    mp_size_t n, nl, nr;

    fmpz_mat22_one(M);

    /* fit everything to xln_len */
    nl = SIZ(xln);
    nr = SIZ(xrn);
    n = FLINT_MAX(nl, nr);

    MPZ_REALLOC(xln, n);
    MPZ_REALLOC(xld, n);
    MPZ_REALLOC(yln, n);
    MPZ_REALLOC(yld, n);
    MPZ_REALLOC(xrn, n);
    MPZ_REALLOC(xrd, n);
    MPZ_REALLOC(yrn, n);
    MPZ_REALLOC(yrd, n);

again:

    xln_len = SIZ(xln);
    xld_len = SIZ(xld);
    xrn_len = SIZ(xrn);
    xrd_len = SIZ(xrd);


    xln_ptr = PTR(xln);
    xld_ptr = PTR(xld);
    xrn_ptr = PTR(xrn);
    xrd_ptr = PTR(xrd);

    yln_ptr = PTR(yln);
    yld_ptr = PTR(yld);
    yrn_ptr = PTR(yrn);
    yrd_ptr = PTR(yrd);

    /* supposed xln > xld > 0 */
    FLINT_ASSERT(xln_len >= xld_len);
    FLINT_ASSERT(xld_len > 0);

    /* supposed xrn > xrd > 0 */
    FLINT_ASSERT(xrn_len >= xrd_len);
    FLINT_ASSERT(xrd_len > 0);

    FLINT_ASSERT(xln_ptr[xln_len-1] != 0);
    FLINT_ASSERT(xld_ptr[xld_len-1] != 0);
    FLINT_ASSERT(xrn_ptr[xrn_len-1] != 0);
    FLINT_ASSERT(xrd_ptr[xrd_len-1] != 0);

    nl = xln_len;
    nr = xrn_len;

    if (nl < 3)
        goto cleanup;

    if (nl != xld_len && nl != xld_len + 1)
        goto cleanup;

    if (nr < 3)
        goto cleanup;

    if (nr != xrd_len && nr != xrd_len + 1)
        goto cleanup;


    if (nl == xld_len + 1)
        xld_ptr[nl-1] = 0;

    if (nr == xrd_len + 1)
        xrd_ptr[nr-1] = 0;


    if ((slong)(xln_ptr[nl-1]) < 0)
    {
        A1 = xln_ptr[nl-1];
        A0 = xln_ptr[nl-2];
        B1 = xld_ptr[nl-1];
        B0 = xld_ptr[nl-2];
    }
    else
    {
        int shift;
        count_leading_zeros(shift, xln_ptr[nl - 1]);
        A1 = MPN_EXTRACT_NUMB(shift, xln_ptr[nl-1], xln_ptr[nl-2]);
        A0 = MPN_EXTRACT_NUMB(shift, xln_ptr[nl-2], xln_ptr[nl-3]);
        B1 = MPN_EXTRACT_NUMB(shift, xld_ptr[nl-1], xld_ptr[nl-2]);
        B0 = MPN_EXTRACT_NUMB(shift, xld_ptr[nl-2], xld_ptr[nl-3]);
    }

    written = _red2(s_temp, A1, A0, B1, B0, m);
    if (written <= 0 || s->length + written > limit)
        goto cleanup;

    if (m->det == 1)
    {
        yln_len = _msub(yln_ptr, m->_22, xln_ptr, m->_12, xld_ptr, nl);
        if (yln_len <= 0)
            goto cleanup;

        yld_len = _msub(yld_ptr, m->_11, xld_ptr, m->_21, xln_ptr, nl);
        if (yld_len <= 0)
            goto cleanup;

        yrn_len = _msub(yrn_ptr, m->_22, xrn_ptr, m->_12, xrd_ptr, nr);
        if (yrn_len <= 0)
            goto cleanup;

        yrd_len = _msub(yrd_ptr, m->_11, xrd_ptr, m->_21, xrn_ptr, nr);
        if (yrd_len <= 0)
            goto cleanup;
    }
    else
    {
        yrn_len = _msub(yrn_ptr, m->_12, xld_ptr, m->_22, xln_ptr, nl);
        if (yrn_len <= 0)
            goto cleanup;

        yrd_len = _msub(yrd_ptr, m->_21, xln_ptr, m->_11, xld_ptr, nl);
        if (yrd_len <= 0)
            goto cleanup;

        yln_len = _msub(yln_ptr, m->_12, xrd_ptr, m->_22, xrn_ptr, nr);
        if (yln_len <= 0)
            goto cleanup;

        yld_len = _msub(yld_ptr, m->_21, xrn_ptr, m->_11, xrd_ptr, nr);
        if (yld_len <= 0)
            goto cleanup;
    }

    /* check yl > 1 */
    if (yln_len <= yld_len)
    {
        if (yln_len != yld_len)
            goto cleanup;

        if (mpn_cmp(yln_ptr, yld_ptr, yln_len) <= 0)
            goto cleanup;
    }

    SIZ(yln) = yln_len;
    SIZ(yld) = yld_len;
    SIZ(yrn) = yrn_len;
    SIZ(yrd) = yrd_len;

    if (needM)
        fmpz_mat22_rmul_ui(M, m);

    fmpz_poly_fit_length(s, s->length + written);
    for (i = 0; i < written; i++)
        fmpz_set_ui(s->coeffs + s->length + i, s_temp[i]);
    s->length += written;
    FLINT_ASSERT(s->length <= limit);

    my_mpz_swap(xln, yln);
    my_mpz_swap(xld, yld);
    my_mpz_swap(xrn, yrn);
    my_mpz_swap(xrd, yrd);
    goto again;

cleanup:

    /* never wrong */
    SIZ(yln) = 0;
    SIZ(yld) = 0;
    SIZ(yrn) = 0;
    SIZ(yrd) = 0;

    return;

}



void fmpq_ball_get_cfrac(
    fmpz_poly_t s, const slong limit,
    fmpz_mat22_t M, int needM,
    fmpq_ball_t x)
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

    FLINT_ASSERT(fmpq_ball_gt_one(x));

again:

/*
    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(x->exact || _fmpq_cmp(x->left_num, x->left_den,
                                       x->right_num, x->right_den) <= 0);
    FLINT_ASSERT(fmpq_ball_gt_one(x));
*/

    if (s->length >= limit)
        goto cleanup;

    k = fmpz_bits(x->left_num);

    if (k > 4*FLINT_BITS)
    {
        if (k > 500*FLINT_BITS)
        {
            k = k/8*5;
            goto split;
        }
        else
        {
            goto lehmer;
        }
    }

gauss:

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    if (s->length >= limit)
        goto cleanup;

    fmpz_fdiv_qr(q, r, x->left_num, x->left_den);
    fmpq_ball_apply_mat22_inv_elem2(y, q, r, x);
    if (!fmpq_ball_gt_one(y))
        goto cleanup;

    fmpq_ball_swap(x, y);

    if (needM)
        fmpz_mat22_rmul_elem(M, q);

    fmpz_poly_fit_length(s, s->length + 1);
    fmpz_swap(s->coeffs + s->length, q);
    s->length++;
    FLINT_ASSERT(s->length <= limit);

    goto gauss;

lehmer:

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    if (s->length >= limit)
        goto cleanup;

    if (x->exact)
    {
        if (COEFF_IS_MPZ(*x->left_num) && COEFF_IS_MPZ(*x->left_den))
        {
            fmpq_ball_get_cfrac_lehmer_exact(s, limit, N, needM,
                       COEFF_TO_PTR(*x->left_num), COEFF_TO_PTR(*x->left_den),
                       _fmpz_promote(y->left_num), _fmpz_promote(y->left_den));
            _fmpz_demote_val(x->left_num);
            _fmpz_demote_val(x->left_den);
            _fmpz_demote_val(y->left_num);
            _fmpz_demote_val(y->left_den);

            if (needM)
                fmpz_mat22_rmul(M, N);
        }
    }
    else
    {
        if (COEFF_IS_MPZ(*x->left_num) && COEFF_IS_MPZ(*x->left_den) &&
            COEFF_IS_MPZ(*x->right_num) && COEFF_IS_MPZ(*x->right_den))
        {
            fmpq_ball_get_cfrac_lehmer_inexact(s, limit, N, needM,
                       COEFF_TO_PTR(*x->left_num), COEFF_TO_PTR(*x->left_den),
                     COEFF_TO_PTR(*x->right_num), COEFF_TO_PTR(*x->right_den),
                       _fmpz_promote(y->left_num), _fmpz_promote(y->left_den),
                     _fmpz_promote(y->right_num), _fmpz_promote(y->right_den));
            _fmpz_demote_val(x->left_num);
            _fmpz_demote_val(x->left_den);
            _fmpz_demote_val(x->right_num);
            _fmpz_demote_val(x->right_den);
            _fmpz_demote_val(y->left_num);
            _fmpz_demote_val(y->left_den);
            _fmpz_demote_val(y->right_num);
            _fmpz_demote_val(y->right_den);

            if (needM)
                fmpz_mat22_rmul(M, N);
        }
    }
/*
    FLINT_ASSERT(x->exact || _fmpq_cmp(x->left_num, x->left_den,
                                       x->right_num, x->right_den) <= 0);
*/

single_step:

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    if (s->length >= limit)
        goto cleanup;

    fmpz_fdiv_qr(q, r, x->left_num, x->left_den);
    fmpq_ball_apply_mat22_inv_elem2(y, q, r, x);
    if (!fmpq_ball_gt_one(y))
        goto cleanup;

    fmpq_ball_swap(x, y);

    if (needM)
        fmpz_mat22_rmul_elem(M, q);

    fmpz_poly_fit_length(s, s->length + 1);
    fmpz_swap(s->coeffs + s->length, q);
    s->length++;
    FLINT_ASSERT(s->length <= limit);

    goto again;

split:

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    fmpq_ball_chop(y, x, k);
    if (!fmpq_ball_gt_one(y))
        goto single_step;

/*    FLINT_ASSERT(fmpq_ball_contains(y, x));*/
    fmpq_ball_get_cfrac(s, limit, N, 1, y);
    if (fmpz_mat22_is_one(N))
        goto single_step;

    fmpq_ball_apply_mat22_inv(y, N, x);

    if (needM)
        fmpz_mat22_rmul(M, N);

    fmpq_ball_swap(x, y);
    fmpq_ball_get_cfrac(s, limit, N, needM, x);

    if (needM)
        fmpz_mat22_rmul(M, N);

cleanup:

    fmpz_clear(q);
    fmpz_clear(r);
    fmpq_ball_clear(y);
    fmpz_mat22_clear(N);

    FLINT_ASSERT(s->length <= limit);
    FLINT_ASSERT(fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

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
    fmpq_ball_t x;
/*
    FLINT_ASSERT(fmpz_sgn(l_den) > 0);
    FLINT_ASSERT(fmpz_sgn(r_den) > 0);
    FLINT_ASSERT(_fmpq_cmp(l_num, l_den, r_num, r_den) <= 0);
*/
    fmpz_init(q);
    fmpz_poly_init(s);
    fmpq_ball_init(x);
    fmpz_mat22_init(M);

    fmpz_mat22_one(M);

    fmpz_set(x->left_num, l_num);
    fmpz_set(x->left_den, l_den);
    fmpz_set(x->right_num, r_num);
    fmpz_set(x->right_den, r_den);

    if (fmpz_cmp(x->left_num, x->left_den) > 0)
    {
        /* 1 < x */
        fmpq_ball_get_cfrac(s, WORD_MAX, M, 1, x);
    }
    else if (fmpz_sgn(x->left_num) > 0
                 && fmpz_cmp(x->right_num, x->right_den) < 0)
    {
        /* 0 < x < 1 */
        fmpz_swap(x->left_den, x->right_num);
        fmpz_swap(x->left_num, x->right_den);
        fmpq_ball_get_cfrac(s, WORD_MAX, M, 1, x);
        fmpz_zero(q);
        fmpz_mat22_lmul_elem(M, q);
    }
    else
    {
        fmpq_ball_t y;
        fmpq_ball_init(y);
        fmpz_fdiv_q(q, x->left_num, x->left_den);
        fmpq_ball_apply_mat22_inv_elem(y, q, x);
        if (fmpq_ball_gt_one(y))
        {
            fmpq_ball_swap(x, y);
            fmpq_ball_get_cfrac(s, WORD_MAX, M, 1, x);
            fmpz_mat22_lmul_elem(M, q);
        }
        fmpq_ball_clear(y);
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
    fmpq_ball_clear(x);
    fmpz_mat22_clear(M);

/*    FLINT_ASSERT(_fmpq_is_canonical(mid_num, mid_den));*/
    return;
}




slong fmpq_get_cfrac(fmpz * c, fmpq_t rem, const fmpq_t f, slong limit)
{
    slong i;
    int cmp, num_sign, den_sign;
    fmpq_ball_t x;
    fmpz_poly_t s;
    fmpz_mat22_t M;
/*
#if WANT_ASSERT
    int input_is_canonical;
#endif
*/
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
/*
#if WANT_ASSERT
    input_is_canonical = fmpq_is_canonical(f);
#endif
*/
    fmpz_mat22_init(M);
    fmpz_mat22_one(M);

    fmpq_ball_init(x);
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
        fmpq_ball_get_cfrac(s, limit, M, 0, x);
    }
    else
    {
        fmpz_poly_fit_length(s, 1);
        if (cmp < 0 && fmpz_sgn(x->left_num) >= 0)
            fmpz_zero(s->coeffs + 0);
        else
            fmpz_fdiv_qr(s->coeffs + s->length, x->left_num, x->left_num, x->left_den);
        s->length = 1;
        fmpz_swap(x->left_num, x->left_den);

        if (!fmpz_is_zero(x->left_den))
            fmpq_ball_get_cfrac(s, limit, M, 0, x);
    }

    FLINT_ASSERT(fmpz_is_zero(x->left_den) || fmpq_ball_gt_one(x));
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
/*
    FLINT_ASSERT(!input_is_canonical || fmpq_is_canonical(rem));
*/
    /* write terms */
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

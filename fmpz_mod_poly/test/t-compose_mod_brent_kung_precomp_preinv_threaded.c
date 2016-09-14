/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>
#include <stdio.h>

#undef ulong

#include <gmp.h>
#include <pthread.h>

#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    
    flint_printf("compose_mod_brent_kung_precomp_preinv_threaded....");
    fflush(stdout);

#if HAVE_PTHREAD && (HAVE_TLS || FLINT_REENTRANT)

    /* check precomputation */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, cinv, * tmp;
        fmpz_t p;
        fmpz_mat_t B, *C;
        slong j, num_threads;
        fmpz_mod_poly_matrix_precompute_arg_t * args1;
        pthread_t *threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_get_num_threads();

        threads = flint_malloc(sizeof(pthread_t) * num_threads);
        tmp = flint_malloc(sizeof(fmpz_mod_poly_t) * num_threads);

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(cinv, p);
        for (j= 0; j < num_threads; j++)
            fmpz_mod_poly_init(tmp[j], p);

        fmpz_mod_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 20) + 1);
        do
        {
            fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);
        } while (c->length < 2);

        fmpz_mod_poly_reverse(cinv, c, c->length);
        fmpz_mod_poly_inv_series_newton(cinv, cinv, c->length);

        fmpz_mat_init(B, n_sqrt(c->length-1)+1, c->length-1);
        fmpz_mod_poly_precompute_matrix(B, b, c, cinv);

        args1 = flint_malloc(sizeof(fmpz_mod_poly_matrix_precompute_arg_t)
                             * num_threads);
        C = flint_malloc(sizeof(fmpz_mat_t) * num_threads);

        for (j = 0; j < num_threads; j++)
        {
            fmpz_mat_init(C[j], n_sqrt(c->length - 1) + 1, c->length - 1);
            fmpz_mod_poly_set(tmp[j], b);
            fmpz_mod_poly_rem(tmp[j], tmp[j], c);
            if (tmp[j]->length < c->length - 1)
            {
                fmpz_mod_poly_fit_length(tmp[j], c->length - 1);
                _fmpz_vec_zero(tmp[j]->coeffs + tmp[j]->length,
                               c->length - 1 - b->length);
            }

            args1[j].A        = *C[j];
            args1[j].poly1    = *tmp[j];
            args1[j].poly2    = *c;
            args1[j].poly2inv = *cinv;

            pthread_create(&threads[j], NULL,
                           _fmpz_mod_poly_precompute_matrix_worker, &args1[j]);
        }
        for (j = 0; j < num_threads; j++)
            pthread_join(threads[j], NULL);

        for (j = 0; j < num_threads; j++)
        {
            if (!fmpz_mat_equal(B, C[j]))
            {
                flint_printf("FAIL (precomputation):\n");
                flint_printf("B:\n"); fmpz_mat_print(B); flint_printf("\n");
                flint_printf("C[j]:\n"); fmpz_mat_print(C[j]); flint_printf("\n");
                flint_printf("a:\n"); fmpz_mod_poly_print(a); flint_printf("\n");
                flint_printf("b:\n"); fmpz_mod_poly_print(b); flint_printf("\n");
                flint_printf("c:\n"); fmpz_mod_poly_print(c); flint_printf("\n");
                abort();
            }
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mat_clear     (B);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(cinv);
        for (j = 0; j < num_threads; j++)
        {
            fmpz_mod_poly_clear(tmp[j]);
            fmpz_mat_clear(C[j]);
        }
        flint_free(C);
        flint_free(tmp);
        flint_free(args1);
        flint_free(threads);
    }

    /* check composition */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, cinv, d, *res;
        fmpz_t p;
        fmpz_mat_t B;
        slong j, num_threads;
        fmpz_mod_poly_compose_mod_precomp_preinv_arg_t * args1;
        pthread_t *threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_get_num_threads();

        threads = flint_malloc(sizeof(pthread_t) * num_threads);
        res = flint_malloc(sizeof(fmpz_mod_poly_t) * num_threads);

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(cinv, p);
        fmpz_mod_poly_init(d, p);
        for (j= 0; j < num_threads; j++)
            fmpz_mod_poly_init(res[j], p);

        fmpz_mod_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1);
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 20) + 1);
        do
        {
            fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);
        } while (c->length < 2);

        fmpz_mod_poly_reverse(cinv, c, c->length);
        fmpz_mod_poly_inv_series_newton(cinv, cinv, c->length);

        fmpz_mat_init(B, n_sqrt(c->length-1)+1, c->length-1);
        fmpz_mod_poly_precompute_matrix(B, b, c, cinv);

        fmpz_mod_poly_rem(a, a, c);
        fmpz_mod_poly_compose_mod(d, a, b, c);

        args1 = flint_malloc(num_threads *
                        sizeof(fmpz_mod_poly_compose_mod_precomp_preinv_arg_t));

        for (j = 0; j < num_threads; j++)
        {
            fmpz_mod_poly_fit_length(res[j], c->length - 1);
            _fmpz_mod_poly_set_length(res[j], c->length - 1);
            args1[j].A        = *B;
            args1[j].res      = *res[j];
            args1[j].poly1    = *a;
            args1[j].poly3    = *c;
            args1[j].poly3inv = *cinv;

            pthread_create(&threads[j], NULL,
                           _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker, &args1[j]);
        }
        for (j = 0; j < num_threads; j++)
        {
            pthread_join(threads[j], NULL);
            _fmpz_mod_poly_normalise(res[j]);
        }

        for (j = 0; j < num_threads; j++)
        {
            if (!fmpz_mod_poly_equal(d, res[j]))
            {
                flint_printf("FAIL (composition):\n");
                flint_printf("res[j]:\n"); fmpz_mod_poly_print(res[j]); flint_printf("\n");
                flint_printf("d:\n"); fmpz_mod_poly_print(d); flint_printf("\n");
                flint_printf("a:\n"); fmpz_mod_poly_print(a); flint_printf("\n");
                flint_printf("b:\n"); fmpz_mod_poly_print(b); flint_printf("\n");
                flint_printf("c:\n"); fmpz_mod_poly_print(c); flint_printf("\n");
                abort();
            }
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mat_clear     (B);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(cinv);
        fmpz_mod_poly_clear(d);
        for (j = 0; j < num_threads; j++)
            fmpz_mod_poly_clear(res[j]);
        flint_free(res);
        flint_free(args1);
        flint_free(threads);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;

#else

   FLINT_TEST_CLEANUP(state);

   flint_printf("SKIPPED\n");
   return 0;

#endif

}


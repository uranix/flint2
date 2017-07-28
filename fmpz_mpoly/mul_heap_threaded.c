/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include <pthread.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"


pthread_mutex_t displock;
typedef struct
{
    slong idx;
    slong nthreads;
    volatile slong stage;
    slong score;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
}
mul_heap_threaded_arg_t;


void * _fmpz_mpoly_mul_heap_threaded_worker(void * arg_ptr)
{
    mul_heap_threaded_arg_t * arg = (mul_heap_threaded_arg_t *) arg_ptr;
    mul_heap_threaded_arg_t *parg = (mul_heap_threaded_arg_t *) arg_ptr - 1;
    pthread_mutex_lock(&displock);
    flint_printf("hello from thread %p idx %wd, prev = %p\n", arg, arg->idx, parg);
    pthread_mutex_unlock(&displock);

    /*
        do first task
        signal when done
    */
    pthread_mutex_lock(&displock);
    flint_printf("thread %wd doing task 1\n", arg->idx);
    pthread_mutex_unlock(&displock);

    arg->stage = 1;
    if (arg->idx + 1 < arg->nthreads)
    {
        pthread_mutex_lock(&arg->mutex);
        pthread_cond_signal(&arg->cond);
        pthread_mutex_unlock(&arg->mutex);
    }

    /*
        wait on previous thread to finish first task
        do second task
    */
    if (arg->idx > 0)
    {
        while (parg->stage < 1)
        {
            pthread_mutex_lock(&parg->mutex);
            pthread_cond_wait(&parg->cond, &parg->mutex);
            pthread_mutex_unlock(&parg->mutex);
        }
    }

    pthread_mutex_lock(&displock);
    flint_printf("thread %wd doing task 2\n", arg->idx);
    pthread_mutex_unlock(&displock);


    pthread_exit(NULL);
}


slong _fmpz_mpoly_mul_heap_threaded(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                 const fmpz * poly3, const ulong * exp3, slong len3,
                                           slong N, ulong maskhi, ulong masklo)
{
    pthread_t * threads;
    mul_heap_threaded_arg_t * args;
    slong i, num_threads;

    num_threads = flint_get_num_threads();
    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(mul_heap_threaded_arg_t) * num_threads);

    pthread_mutex_init(&displock, NULL);

    for (i = 0; i < num_threads; i++)
    {
        args[i].idx = i;
        args[i].nthreads = num_threads;
        args[i].stage = 0;
        pthread_mutex_init(&args[i].mutex, NULL);
        pthread_cond_init(&args[i].cond, NULL);

        pthread_create(&threads[i], NULL, _fmpz_mpoly_mul_heap_threaded_worker, &args[i]);
    }

    for (i = 0; i < num_threads; i++)
        pthread_join(threads[i], NULL);

    flint_free(threads);
    flint_free(args);

    return 0;
}

void fmpz_mpoly_mul_heap_threaded(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
   slong i, bits, exp_bits, N, len = 0;
   ulong * max_degs2;
   ulong * max_degs3;
   ulong maskhi, masklo;
   ulong max;
   ulong * exp2 = poly2->exps, * exp3 = poly3->exps;
   int free2 = 0, free3 = 0;

   TMP_INIT;

   /* one of the input polynomials is zero */
   if (poly2->length == 0 || poly3->length == 0)
   {
      fmpz_mpoly_zero(poly1, ctx);

      return;
   }

   TMP_START;

   /* compute maximum degree of any variable */
   max_degs2 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));
   max_degs3 = (ulong *) TMP_ALLOC(ctx->n*sizeof(ulong));

   fmpz_mpoly_max_degrees(max_degs2, poly2, ctx);
   fmpz_mpoly_max_degrees(max_degs3, poly3, ctx);

   max = 0;

   for (i = 0; i < ctx->n; i++)
   {
      max_degs3[i] += max_degs2[i];
      /*check exponents won't overflow */
      if (max_degs3[i] < max_degs2[i] || 0 > (slong) max_degs3[i]) 
         flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_mul_johnson");

      if (max_degs3[i] > max)
         max = max_degs3[i];
   }

   /* compute number of bits to store maximum degree */
   bits = FLINT_BIT_COUNT(max);
   if (bits >= FLINT_BITS)
      flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_mul_johnson");

   exp_bits = 8;
   while (bits >= exp_bits) /* extra bit required for signs */
      exp_bits *= 2;

   exp_bits = FLINT_MAX(exp_bits, poly2->bits);
   exp_bits = FLINT_MAX(exp_bits, poly3->bits);

   masks_from_bits_ord(maskhi, masklo, exp_bits, ctx->ord);
   /* number of words exponent vectors packed into */
   N = (exp_bits*ctx->n - 1)/FLINT_BITS + 1;

   /* ensure input exponents are packed into same sized fields as output */
   if (exp_bits > poly2->bits)
   {
      free2 = 1;
      exp2 = (ulong *) flint_malloc(N*poly2->length*sizeof(ulong));
      mpoly_unpack_monomials(exp2, exp_bits, poly2->exps, poly2->bits,
                                                        poly2->length, ctx->n);
   }

   if (exp_bits > poly3->bits)
   {
      free3 = 1;
      exp3 = (ulong *) flint_malloc(N*poly3->length*sizeof(ulong));
      mpoly_unpack_monomials(exp3, exp_bits, poly3->exps, poly3->bits,
                                                        poly3->length, ctx->n);
   }

   /* deal with aliasing and do multiplication */
   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(temp, exp_bits, ctx);
      temp->bits = exp_bits;

      /* algorithm more efficient if smaller poly first */
      if (poly2->length >= poly3->length)
         len = _fmpz_mpoly_mul_heap_threaded(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                               N, maskhi, masklo);
      else
         len = _fmpz_mpoly_mul_heap_threaded(&temp->coeffs, &temp->exps, &temp->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                               N, maskhi, masklo);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length - 1, ctx);
      fmpz_mpoly_fit_bits(poly1, exp_bits, ctx);
      poly1->bits = exp_bits;

      /* algorithm more efficient if smaller poly first */
      if (poly2->length > poly3->length)
         len = _fmpz_mpoly_mul_heap_threaded(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                               N, maskhi, masklo);
      else
         len = _fmpz_mpoly_mul_heap_threaded(&poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly2->coeffs, exp2, poly2->length,
                                      poly3->coeffs, exp3, poly3->length,
                                               N, maskhi, masklo);
   }

   if (free2)
      flint_free(exp2);

   if (free3)
      flint_free(exp3);

   _fmpz_mpoly_set_length(poly1, len, ctx);

   TMP_END;
}

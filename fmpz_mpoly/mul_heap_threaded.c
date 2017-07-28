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
    ulong * dexp;
    slong * dexp_ind;
    slong dexp_score;
    slong score;
    fmpz * coeff2; ulong * exp2; slong len2;
    fmpz * coeff3; ulong * exp3; slong len3;
    slong len1;
    slong alloc1;
    ulong * exp1;
    fmpz * coeff1;

    slong N;
    ulong maskhi; ulong masklo;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
}
mul_heap_threaded_arg_t;



/*
    Set poly1 to poly2*poly3 using Johnson's heap method. The function
    realocates its output and returns the length of the product. This
    version of the function assumes the exponent vectors take N words.
    only terms starting at indicies in "start" and strictly greater than
    endexp are written
*/
slong _fmpz_mpoly_mul_heap_part(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                 const fmpz * poly3, const ulong * exp3, slong len3,
            slong * start, slong * end, slong N, ulong maskhi, ulong masklo)
{
    slong i, k;
    slong next_free, Q_len = 0, heap_len = 2; /* heap zero index unused */
    mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    mpoly_heap_t ** Q;
    mpoly_heap_t * x;
    fmpz * p1 = *poly1;
    ulong * e1 = *exp1;
    ulong cy;
    ulong c[3], p[2]; /* for accumulating coefficients */
    ulong * exp, * exps;
    ulong ** exp_list;
    slong exp_next;
    int first, small;

   TMP_INIT;

/*    pthread_mutex_lock(&displock);*/
    for (i = 0; i < len2; i++)
    {
        flint_printf("a index %wd: start %wd end %wd\n",i,start[i],end[i]);
    }
/*    pthread_mutex_unlock(&displock);*/

    return 0;



   TMP_START;

   /* whether input coeffs are small, thus output coeffs fit in three words */
   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

   heap = (mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap_s));
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (mpoly_heap_t **) TMP_ALLOC(len2*sizeof(mpoly_heap_t *));
   /* allocate space for exponent vectors of N words */
   exps = (ulong *) TMP_ALLOC(len2*N*sizeof(ulong));
   /* list of pointers to allocated exponent vectors */
   exp_list = (ulong **) TMP_ALLOC(len2*sizeof(ulong *));

   for (i = 0; i < len2; i++)
      exp_list[i] = exps + i*N;

   /* start with no heap nodes and no exponent vectors in use */
   next_free = 0;
   exp_next = 0;

   /* put (0, 0, exp2[0] + exp3[0]) on heap */
   x = chain + next_free++;
   x->i = 0;
   x->j = 0;
   x->next = NULL;

   heap[1].next = x;
   heap[1].exp = exp_list[exp_next++];

   mpoly_monomial_add(heap[1].exp, exp2, exp3, N);




         /* take node from store */
/*
         x = Q[--Q_len];
     
         if (x->j == 0 && x->i < len2 - 1)
         {
            mpoly_heap_t * x2 = chain + next_free++;
            x2->i = x->i + 1;
            x2->j = 0;
            x2->next = NULL;
*/
            /* insert (x->i + 1, 0, exps[x->i + 1] + exp3[0]) in heap */
/*            
            mpoly_monomial_add(exp_list[exp_next], exp2 + (x->i + 1)*N, exp3, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x2, &heap_len,
                                                            N, maskhi, masklo))
               exp_next--;
         }
*/


   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);

   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get pointer to exponent field of heap top */
      exp = heap[1].exp;
      
      /* realloc output poly ready for next product term */
      k++;
      _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, N);

      /* whether we are on first coeff product for this output exponent */
      first = 1;

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && mpoly_monomial_equal(heap[1].exp, exp, N))
      {
         /* pop chain from heap and set exponent field to be reused */
         exp_list[--exp_next] = heap[1].exp;

         x = _mpoly_heap_pop(heap, &heap_len, N, maskhi, masklo);
         
         /* if output coeffs will fit in three words */
         if (small)
         {
            /* compute product of input poly coeffs */
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));

               /* set output monomial */
               mpoly_monomial_set(e1 + k*N, exp, N);

               first = 0; 
            } else /* addmul product of input poly coeffs */
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
            }
      
            /* temporarily store pointer to this node */
            if (x->j < len3 - 1 || x->j == 0)
               Q[Q_len++] = x;

            /* for every node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               /* temporarily store pointer to this node */
               if (x->j < len3 - 1 || x->j == 0)
                  Q[Q_len++] = x;
            }
         } else /* output coeffs require multiprecision */
         {
            if (first) /* compute product of input poly coeffs */
            {
               fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
               /* set output monomial */
               mpoly_monomial_set(e1 + k*N, exp, N);

               first = 0; 
            } else /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
      
            /* temporarily store pointer to this node */
            if (x->j < len3 - 1 || x->j == 0)
               Q[Q_len++] = x;

            /* for each node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

               /* temporarily store pointer to this node */
               if (x->j < len3 - 1 || x->j == 0)
                  Q[Q_len++] = x;
            }
         }
      }
      
      /* for each node temporarily stored */
      while (Q_len > 0)
      {
         /* take node from store */
         x = Q[--Q_len];
     
         if (x->j == 0 && x->i < len2 - 1)
         {
            mpoly_heap_t * x2 = chain + next_free++;
            x2->i = x->i + 1;
            x2->j = 0;
            x2->next = NULL;

            /* insert (x->i + 1, 0, exps[x->i + 1] + exp3[0]) in heap */
            
            mpoly_monomial_add(exp_list[exp_next], exp2 + (x->i + 1)*N, exp3, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x2, &heap_len,
                                                            N, maskhi, masklo))
               exp_next--;
         }

         if (x->j < len3 - 1)
         {
            x->j++;
            x->next = NULL;

            /* insert (x->i, x->j + 1, exps[x->i] + exp3[x->j + 1]) in heap */

            mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N, exp3 + x->j*N, N);

            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len,
                                                            N, maskhi, masklo))
               exp_next--;
         }
      }     
   
      /* set output poly coeff from temporary accumulation, if not multiprec */
      if (small)
         fmpz_set_signed_uiuiui(p1 + k, c[2], c[1], c[0]);

      if (fmpz_is_zero(p1 + k))
         k--;
   }

   k++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}









void * _fmpz_mpoly_mul_heap_threaded_worker(void * arg_ptr)
{
    mul_heap_threaded_arg_t * arg = (mul_heap_threaded_arg_t *) arg_ptr;
    mul_heap_threaded_arg_t *parg = (mul_heap_threaded_arg_t *) arg_ptr - 1;
    slong i, lower, upper, * dummy;

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

    lower = ((arg->idx + 1)*arg->len2*arg->len3) / arg->nthreads;
    upper = lower;

    fmpz_mpoly_search_monomials(arg->dexp, &arg->dexp_score, arg->dexp_ind,
                            lower, upper,
                            arg->exp2, arg->len2, arg->exp3, arg->len3,
                                     arg->N, arg->maskhi, arg->masklo);

    pthread_mutex_lock(&displock);
    flint_printf("thread %wd finished task 1\n", arg->idx);
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


    if (arg->idx == 0)
    {
        dummy = flint_malloc(arg->len2*sizeof(slong));
        for (i = 0; i < arg->len2; i++)
            dummy[i] = arg->len3;

    }
    arg->len1 = _fmpz_mpoly_mul_heap_part(
                     &arg->coeff1, &arg->exp1, &arg->alloc1,
                      arg->coeff2,  arg->exp2,  arg->len2,
                      arg->coeff3,  arg->exp3,  arg->len3,
                       arg->dexp_ind, ((arg->idx == 0) ? dummy : (parg->dexp_ind)),
                                             arg->N, arg->maskhi, arg->masklo);
    if (arg->idx == 0)
        flint_free(dummy);

    pthread_mutex_unlock(&displock);

    pthread_exit(NULL);
}


slong _fmpz_mpoly_mul_heap_threaded(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * coeff2, const ulong * exp2, slong len2,
                 const fmpz * coeff3, const ulong * exp3, slong len3,
                                           slong N, ulong maskhi, ulong masklo)
{
    pthread_t * threads;
    mul_heap_threaded_arg_t * args;
    slong i, num_threads;

    num_threads = flint_get_num_threads();
    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(mul_heap_threaded_arg_t) * num_threads);

    pthread_mutex_init(&displock, NULL);

flint_printf("len2: %wd  len3: %wd\n",len2,len3);

    for (i = 0; i < num_threads; i++)
    {
        args[i].idx = i;
        args[i].nthreads = num_threads;
        args[i].stage = 0;
        args[i].dexp = (ulong *) flint_malloc(N*sizeof(ulong));
        args[i].dexp_ind = (slong *) flint_malloc(len2*sizeof(slong));

        args[i].coeff2 = coeff2;
        args[i].exp2 = exp2;
        args[i].len2 = len2;
        args[i].coeff3 = coeff3;
        args[i].exp3 = exp3;
        args[i].len3 = len3;

        args[i].N = N;
        args[i].maskhi = maskhi;
        args[i].masklo = masklo;
        args[i].len1 = 0;
        args[i].alloc1 = len2 + len3/num_threads;
        args[i].exp1 = (ulong *) flint_malloc(args[i].alloc1*N*sizeof(ulong)); 
        args[i].coeff1 = (fmpz *) flint_calloc(args[i].alloc1, sizeof(fmpz));

        pthread_mutex_init(&args[i].mutex, NULL);
        pthread_cond_init(&args[i].cond, NULL);

        pthread_create(&threads[i], NULL, _fmpz_mpoly_mul_heap_threaded_worker, &args[i]);
    }

    for (i = num_threads - 1; i >= 0; i--)
    {
        pthread_join(threads[i], NULL);
    }


    for (i = num_threads - 1; i >= 0; i--)
    {
        flint_free(args[i].coeff1);
        flint_free(args[i].exp1);
        flint_free(args[i].dexp_ind);
        flint_free(args[i].dexp);
    }

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

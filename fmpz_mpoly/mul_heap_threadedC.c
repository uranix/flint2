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
#include <assert.h>
#include "profiler.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"


/*
    Set poly1 to poly2*poly3 using Johnson's heap method. The function
    realocates its output and returns the length of the product. This
    version of the function assumes the exponent vectors all fit in a
    single word. Assumes input polys are nonzero.
    Only terms t with start >= t > end are written;
    "start" and "end" are not monomials but arrays of indicides into exp3
*/
slong _fmpz_mpoly_mul_heap_part1C(fmpz ** poly1, ulong ** exp1, slong * alloc,
              const fmpz * poly2, const ulong * exp2, slong len2,
              const fmpz * poly3, const ulong * exp3, slong len3,
                                      slong * start, slong * end, ulong maskhi)
{
    slong i, k;
    slong next_free, Q_len = 0, heap_len = 1; /* heap zero index unused */
    mpoly_heap1_s * heap;
    mpoly_heap_t * chain;
    mpoly_heap_t ** Q;
    mpoly_heap_t * x;
    fmpz * p1 = *poly1;
    ulong * e1 = *exp1;
    ulong exp, cy;
    ulong c[3], p[2]; /* for accumulating coefficients */
    int first, small;
    TMP_INIT;

    TMP_START;

    /* whether input coeffs are small, thus output coeffs fit in three words */
    small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

    heap = (mpoly_heap1_s *) TMP_ALLOC((len2 + 1)*sizeof(mpoly_heap1_s));
    /* alloc array of heap nodes which can be chained together */
    chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
    /* space for temporary storage of pointers to heap nodes */
    Q = (mpoly_heap_t **) TMP_ALLOC(len2*sizeof(mpoly_heap_t *));
   
    /* start with no heap nodes in use */
    next_free = 0;

    /* put all the starting nodes on the heap */
    for (i = 0; i < len2; i++)
        if (start[i] < end[i])
        {
            x = chain + next_free++;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x, &heap_len,
                                                                       maskhi);
        }

    /* output poly index starts at -1, will be immediately updated to 0 */
    k = -WORD(1);

    /* while heap is nonempty */
    while (heap_len > 1)
    {
        /* get exponent field of heap top */
        exp = heap[1].exp;
      
        /* realloc output poly ready for next product term */
        k++;
        _fmpz_mpoly_fit_length(&p1, &e1, alloc, k + 1, 1);

        /* whether we are on first coeff product for this output exponent */
        first = 1;

        /* set temporary coeff to zero */
        c[0] = c[1] = c[2] = 0;

        /* while heap nonempty and contains chain with current output exponent */
        while (heap_len > 1 && heap[1].exp == exp)
        {
            /* pop chain from heap */
            x = _mpoly_heap_pop1(heap, &heap_len, maskhi);

            /* if output coeffs will fit in three words */
            if (small)
            {
                /* compute product of input poly coeffs */
                if (first)
                {
                    smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
                    c[2] = -(c[1] >> (FLINT_BITS - 1));

                    /* set output monomial */
                    e1[k] = exp;
                    first = 0; 
                } else /* addmul product of input poly coeffs */
                {
                    smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
                    add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
                    c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
                }
      
                /* temporarily store pointer to this node */
                if (x->j < len3 - 1)
                    Q[Q_len++] = x;

                /* for every node in this chain */
                while ((x = x->next) != NULL)
                {          
                    /* addmul product of input poly coeffs */
                    smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
                    add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
                    c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

                    /* temporarily store pointer to this node */
                    if (x->j < len3 - 1)
                        Q[Q_len++] = x;
                }
           } else /* output coeffs require multiprecision */
           {
                if (first) /* compute product of input poly coeffs */
                {
                    fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
                    e1[k] = exp;
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

            if (x->j + 1 < end[x->i])
            {
                x->j++;
                x->next = NULL;
                _mpoly_heap_insert1(heap, exp2[x->i] + exp3[x->j], x, &heap_len,
                                                                       maskhi);
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



/*
    Set poly1 to poly2*poly3 using Johnson's heap method. The function
    realocates its output and returns the length of the product. This
    version of the function assumes the exponent vectors take N words.
    Only terms t with start >= t > end are written;
    "start" and "end" are not monomials but arrays of indicides into exp3
*/
slong _fmpz_mpoly_mul_heap_partC(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * poly2, const ulong * exp2, slong len2,
                 const fmpz * poly3, const ulong * exp3, slong len3,
            slong * start, slong * end, slong N, ulong maskhi, ulong masklo)
{
    slong i, k;
    slong next_free, Q_len = 0, heap_len = 1; /* heap zero index unused */
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

    /* if exponent vectors fit in single word, call special version */
    if (N == 1)
        return _fmpz_mpoly_mul_heap_part1C(poly1, exp1, alloc,
                                          poly2, exp2, len2,
                                          poly3, exp3, len3,
                                                    start, end, maskhi);

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
    {
      exp_list[i] = exps + i*N;
    }

    /* start with no heap nodes and no exponent vectors in use */
    next_free = 0;
    exp_next = 0;


    /* put all the starting nodes on the heap */
    for (i = 0; i < len2; i++)
        if (start[i] < end[i])
        {
            x = chain + next_free++;
            x->i = i;
            x->j = start[i];
            x->next = NULL;
            mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N,
                                                             exp3 + x->j*N, N);
            if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x, &heap_len,
                                                            N, maskhi, masklo))
               exp_next--;
        }

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

            /* do not put next on heap if it is past the end */
            if (x->j + 1 < end[x->i])
            {
                x->j++;
                x->next = NULL;

                mpoly_monomial_add(exp_list[exp_next], exp2 + x->i*N,
                                                       exp3 + x->j*N, N);
                if (!_mpoly_heap_insert(heap, exp_list[exp_next++], x,
                                                 &heap_len, N, maskhi, masklo))
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


/*
    The workers receive their data via mul_heap_threaded_arg_t.
    The i^th worker (i = 0, ... , n - 1) computes a monomial that divides
    the output into roughly equal portions. This computation is trivial for
    i = n - 1 as it simply grabs the product of the leading monomials in poly2
    and poly . The i^th worker then computes terms that lie between the i^th
    monomial and the (i-1)^st monomial. The i = 0 worker does not need any
    computation to know where to stop and hence does not need to wait on any
    other workers.
*/

typedef struct
{
    pthread_mutex_t mutex;
    slong nthreads;
    slong adivs, bdivs;
    const fmpz * coeff2; const ulong * exp2; slong len2;
    const fmpz * coeff3; const ulong * exp3; slong len3;
    slong N;
    ulong maskhi; ulong masklo;    
    volatile int idx;
}
mul_heap_threadedC_base_t;

typedef struct
{
    slong eat;
    slong eat_idx;
    slong apstart, aplen;
    slong * t1, * t2, *t3;
    slong * line;
    ulong * exp;
    slong score;
    slong lower;
    slong upper;
    slong threadidx;
    slong time;
    slong len1;
    slong alloc1;
    ulong * exp1;
    fmpz * coeff1;
}
mul_heap_threadedC_div_t;


typedef struct
{
    slong idx;
    slong time;
    mul_heap_threadedC_base_t * basep;
    mul_heap_threadedC_div_t * divp;
    volatile slong stage;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
}
mul_heap_threadedC_arg_t;




void * _fmpz_mpoly_mul_heap_threadedC_worker(void * arg_ptr)
{
    mul_heap_threadedC_arg_t * arg = (mul_heap_threadedC_arg_t *) arg_ptr;
    mul_heap_threadedC_div_t * divs;
    mul_heap_threadedC_base_t * base;
    slong k, m, * dummy;
/*
    timeit_t time, time2;
*/

    divs = arg->divp;
    base = arg->basep;
/*
    timeit_start(time);
*/
    pthread_mutex_lock(&base->mutex);
    k = base->idx - 1;
    base->idx = k;
    pthread_mutex_unlock(&base->mutex);

    while (k >= 0)
    {
/*
    pthread_mutex_lock(&base->mutex);
*/
/*
        divs[i].threadidx = arg->idx;
        timeit_start(time2);
*/

        if (k >= base->adivs)
        {
            divs[k].len1 = _fmpz_mpoly_mul_heap_partC(
                         &divs[k].coeff1, &divs[k].exp1, &divs[k].alloc1,
                          base->coeff2 + divs[k].apstart,  base->exp2 + base->N * divs[k].apstart,  divs[k].aplen,
                          base->coeff3,  base->exp3,  base->len3,
                           divs[k].line, divs[k - base->adivs].line,
                                                 base->N, base->maskhi, base->masklo);
        } else
        {
            dummy = flint_malloc(base->len2*sizeof(slong));
            for (m = 0; m < base->len2; m++)
                dummy[m] = base->len3;

            divs[k].len1 = _fmpz_mpoly_mul_heap_partC(
                         &divs[k].coeff1, &divs[k].exp1, &divs[k].alloc1,
                          base->coeff2 + divs[k].apstart,  base->exp2 + base->N * divs[k].apstart, divs[k].aplen,
                          base->coeff3,  base->exp3,  base->len3,
                           divs[k].line, dummy,
                                                 base->N, base->maskhi, base->masklo);
            flint_free(dummy);
        }
/*
        timeit_stop(time2);
        divs[i].time = time2->wall;
*/
/*
    pthread_mutex_unlock(&base->mutex);
*/

        pthread_mutex_lock(&base->mutex);
        k = base->idx - 1;
        base->idx = k;
        pthread_mutex_unlock(&base->mutex);
    }
/*
    timeit_stop(time);
    arg->time = time->wall;
*/
    flint_cleanup();
    return NULL;
}






/*
    start the workers from low to high and join from high to low
    after joining a worker, the output is streamed to poly1
*/
slong _fmpz_mpoly_mul_heap_threadedC(fmpz ** poly1, ulong ** exp1, slong * alloc,
                 const fmpz * coeff2, const ulong * exp2, slong len2,
                 const fmpz * coeff3, const ulong * exp3, slong len3,
                                           slong N, ulong maskhi, ulong masklo)
{
    slong alen = len2, blen = len3;
    slong i, j, k, l, m;
    pthread_t * threads;
    mul_heap_threadedC_arg_t * args;
    mul_heap_threadedC_base_t * base;
    mul_heap_threadedC_div_t * divs;
    fmpz * p1 = * poly1;
    ulong * e1 = * exp1;
    ulong * max_exp;
    int first, not_done;

/*
    timeit_t time;
*/
    base = flint_malloc(sizeof(mul_heap_threadedC_base_t));
    base->nthreads = flint_get_num_threads();
    base->bdivs    = base->nthreads*2;
    base->adivs    = alen < 4 ? 1 : alen < 12 ? 2 : 3;
    base->coeff2 = coeff2;
    base->exp2 = exp2;
    base->len2 = len2;
    base->coeff3 = coeff3;
    base->exp3 = exp3;
    base->len3 = len3;
    base->N = N;
    base->maskhi = maskhi;
    base->masklo = masklo;
    base->idx = base->adivs * base->bdivs;

    divs    = flint_malloc(sizeof(mul_heap_threadedC_div_t) * base->adivs * base->bdivs);
    threads = flint_malloc(sizeof(pthread_t) * base->nthreads);
    args    = flint_malloc(sizeof(mul_heap_threadedC_arg_t) * base->nthreads);

    for (i = 0; i < base->bdivs; i++)
    {
    for (j = 0; j < base->adivs; j++)
    {
/*flint_printf("i: %wd  j: %wd\n",i,j);*/

        k = j + i * base->adivs;

        divs[k].apstart = j * alen / base->adivs;
        divs[k].aplen = (j + 1) * alen / base->adivs - divs[k].apstart;

        divs[k].lower = ((i + 1)*(i + 1)*divs[k].aplen*blen) / (base->bdivs*base->bdivs);
        divs[k].upper = divs[k].lower;

        divs[k].line = NULL;
        divs[k].exp = (ulong *) flint_malloc(N*sizeof(ulong));
        divs[k].t1 = (slong *) flint_malloc(alen*sizeof(slong));
        divs[k].t2 = (slong *) flint_malloc(alen*sizeof(slong));
        divs[k].t3 = (slong *) flint_malloc(alen*sizeof(slong));
        divs[k].eat = 0;

        divs[k].len1   = 0;
        /* everything writes to a new worker poly */
        divs[k].alloc1 = 1 + alen/base->adivs + blen/base->bdivs;
        divs[k].exp1 = (ulong *) flint_malloc(divs[i].alloc1*N*sizeof(ulong)); 
        divs[k].coeff1 = (fmpz *) flint_calloc(divs[i].alloc1, sizeof(fmpz));
    }
    }

/*flint_printf("done alloc divs\n");*/

    for (k = 0; k < base->adivs * base->bdivs; k++)
    {
/*
        timeit_start(time);
*/
        mpoly_search_monomialsB(
            &divs[k].line, divs[k].exp, &divs[k].score,
                            divs[k].t1, divs[k].t2, divs[k].t3,
                            divs[k].lower, divs[k].upper,
                            base->exp2 + N*divs[k].apstart, divs[k].aplen, base->exp3, base->len3,
                                     base->N, base->maskhi, base->masklo);
/*
        timeit_stop(time);
        flint_printf("search %wd: %p score %wd, time %wd\n", i, divs[i].line, divs[i].score, time->wall);
*/
/*        flint_printf("search %wd: %p score %wd\n", k, divs[k].line, divs[k].score);*/

    }




    pthread_mutex_init(&base->mutex, NULL);
    for (i = 0; i < base->nthreads; i++)
    {
        args[i].idx = i;
        args[i].basep = base;
        args[i].divp = divs;
        pthread_create(&threads[i], NULL, _fmpz_mpoly_mul_heap_threadedC_worker, &args[i]);
    }
    for (i = 0; i < base->nthreads; i++)
    {
        pthread_join(threads[i], NULL);
/*
        flint_printf("thread %wd time %wd\n",i,args[i].time);
*/
    }
    pthread_mutex_destroy(&base->mutex);



/*
for (k = base->adivs * base->bdivs - 1; k >= 0; k--)
{
for (i = 0; i < divs[k].len1; i++)
{
    flint_printf("div %wd coeff %wx, exponent %wx\n", k, *((ulong*)(divs[k].coeff1 + i)),*(divs[k].exp1 + N*i));
}
}
*/

    k = 0;
    max_exp = (ulong *) flint_malloc(N*sizeof(ulong));

    /* initialize tops */
    for (j = base->adivs - 1; j >= 0; j--)
    {
        /* find j with largest term on top */
        divs[j].eat_idx = j + (base->bdivs - 1) * base->adivs;
        /* all divs[-].eat have already been set to zero */
        l = divs[j].eat_idx;
        while (l >= 0 && divs[l].eat >= divs[l].len1)
        {
            for (m = divs[l].eat; m < divs[l].alloc1; m++)
                fmpz_clear(divs[l].coeff1 + m);
            flint_free(divs[l].coeff1);
            flint_free(divs[l].exp1);
            flint_free(divs[l].t3);
            flint_free(divs[l].t2);
            flint_free(divs[l].t1);
            flint_free(divs[l].exp);
            l = l - base->adivs;
            divs[j].eat_idx = l;
            /* divs[l].eat is zero because it was initialized to zero */
        }
    }

    do {
        not_done = 0;
        for (m = 0; m < N; m++)
            max_exp[m] = 0;
        for (j = base->adivs - 1; j >= 0; j--)
        {
            /* find j with largest term on top */
            i = divs[j].eat_idx;
            if (i >= 0 && divs[i].eat < divs[i].len1)
            {
                not_done = 1;
                if (mpoly_monomial_lt(divs[i].exp1 + N*divs[i].eat, max_exp, N, maskhi, masklo))
                    mpoly_monomial_set(max_exp, divs[i].exp1 + N*divs[i].eat, N);
            }
        }



        if (not_done)
        {
            first = 1;
            for (j = base->adivs - 1; j >= 0; j--)
            {
                /* find j with largest term on top */
                i = divs[j].eat_idx;
                if (i >= 0 && divs[i].eat < divs[i].len1)
                {
                    if (mpoly_monomial_equal(max_exp, divs[i].exp1 + N*divs[i].eat, N))
                    {
                        _fmpz_mpoly_fit_length(&p1, &e1, alloc, k+1, N);
                        if (first)
                        {
                            fmpz_swap(p1 + k, divs[i].coeff1 + divs[i].eat);
                            mpoly_monomial_set(e1 + N*k, divs[i].exp1 + N*divs[i].eat, N);
                            first = 0;
                        } else {
                            fmpz_add(p1 + k, p1 + k, divs[i].coeff1 + divs[i].eat);
                        }
                        fmpz_clear(divs[i].coeff1 + divs[i].eat);
                        divs[i].eat++;

                        /* all divs[-].eat have already been set to zero */
                        l = divs[j].eat_idx;
                        while (l >= 0 && divs[l].eat >= divs[l].len1)
                        {
                            for (m = divs[l].eat; m < divs[l].alloc1; m++)
                                fmpz_clear(divs[l].coeff1 + m);
                            flint_free(divs[l].coeff1);
                            flint_free(divs[l].exp1);
                            flint_free(divs[l].t3);
                            flint_free(divs[l].t2);
                            flint_free(divs[l].t1);
                            flint_free(divs[l].exp);
                            l = l - base->adivs;
                            divs[j].eat_idx = l;
                            /* divs[l].eat is zero because it was initialized to zero */
                        }
                    }
                }
            }
            assert(first == 0);
            if (fmpz_is_zero(p1 + k))
                fmpz_clear(p1 + k);
            else
                k++;

        }
    } while (not_done);

    flint_free(max_exp);

    flint_free(args);
    flint_free(threads);
    flint_free(divs);
    flint_free(base);

    *poly1 = p1;
    *exp1  = e1;
    return k;
}

void fmpz_mpoly_mul_heap_threadedC(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
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
            len = _fmpz_mpoly_mul_heap_threadedC(
                                    &temp->coeffs, &temp->exps, &temp->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                               N, maskhi, masklo);
        else
            len = _fmpz_mpoly_mul_heap_threadedC(
                                   &temp->coeffs, &temp->exps, &temp->alloc,
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
            len = _fmpz_mpoly_mul_heap_threadedC(
                                &poly1->coeffs, &poly1->exps, &poly1->alloc,
                                      poly3->coeffs, exp3, poly3->length,
                                      poly2->coeffs, exp2, poly2->length,
                                               N, maskhi, masklo);
        else
            len = _fmpz_mpoly_mul_heap_threadedC(
                                &poly1->coeffs, &poly1->exps, &poly1->alloc,
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

/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "flint.h"
#include "ulong_extras.h"
#include "hashmap.h"

typedef struct apow {
    ulong k;
    ulong ak;
} apow_t;

typedef struct {
    ulong n;
    double ninv;
    ulong m;
    ulong am;
    apow_t * table;
} bsgs_struct;

typedef bsgs_struct bsgs_t[1];

static int
apow_cmp(const apow_t * x, const apow_t * y)
{
    return (x->ak < y->ak) ? -1 : (x->ak > y->ak);
}

/* set size of table m=sqrt(nk) to compute k logs in a group of size n */
void
bsgs_table_init(bsgs_t t, ulong a, ulong n, ulong m)
{
    ulong k, ak;
    t->table = (apow_t *)flint_malloc(m * sizeof(apow_t));

    t->n = n;
    t->ninv = n_precompute_inverse(n);
    t->m = m;

    for (k = 0, ak = 1; k < m; k++)
    {
        t->table[k].k = k;
        t->table[k].ak = ak;
        ak = n_mulmod_precomp(ak, a, n, t->ninv);
    }

    t->am = n_invmod(ak, n);
    qsort(t->table, m, sizeof(apow_t), (int(*)(const void*,const void*))apow_cmp);
}

void
bsgs_table_clear(bsgs_t t)
{
    flint_free(t->table);
}

ulong
n_discrete_log_bsgs_table(const bsgs_t t, ulong b)
{
    ulong i;
    apow_t c, * x;

    c.k = 0;
    c.ak = b;
    for (i = 0; i < t->m; i++)
    {
        x = bsearch(&c, t->table, t->m, sizeof(apow_t),
            (int(*)(const void*,const void*))apow_cmp);
        if (x != NULL)
            return i * t->m + x->k;
        c.ak = n_mulmod_precomp(c.ak, t->am, t->n, t->ninv);
    }
    flint_printf("Exception (n_discrete_log_bsgs).  discrete log not found.\n");
    flint_abort();
    return 0; /* not reached, but silence the compiler */
}

ulong
n_discrete_log_bsgs(ulong b, ulong a, ulong n)
{
    ulong m;
    bsgs_t table;

    m = ceil(sqrt((double) n));
    bsgs_table_init(table, a, n, m);
    m = n_discrete_log_bsgs_table(table, b);
    bsgs_table_clear(table);

    return m;
}

/*
   b = array of values to take logs of
   a = base of logarithms
   num = number of values
   ord = order of the group, or 0 if n - 1 should be used
   n = prime
   r = result : array of logarithms, i.e.
   r[i] log_a(b[i]) mod p
*/
void
n_discrete_log_multi_bsgs(ulong * r, const ulong * b, ulong a,
                                                  ulong num, ulong ord, ulong n)
{
    slong i, j;
    ulong num_gs, num_bs, gs_base, gs_pow, gs_shoup;
    ulong * bs, * bs_shoup;
    ulong val, key;
    hashmap1_t gs;

    if (ord == 0)
       ord = n - 1;

    /* number of giant steps */
    num_gs = ceil(sqrt((double) ord))*ceil(sqrt((double) num));

    /* number of baby steps */
    num_bs = (ord - 1)/num_gs + 1;

    bs = flint_malloc(num_bs*sizeof(ulong));
    bs_shoup = flint_malloc(num_bs*sizeof(ulong));

    bs[0] = a;
    bs_shoup[0] = n_mulmod_precomp_shoup(a, n);
    for (i = 1; i < num_bs; i++)
    {
        bs[i] = n_mulmod_shoup(a, bs[i - 1], bs_shoup[0], n);
        bs_shoup[i] = n_mulmod_precomp_shoup(bs_shoup[i], n);
    }

    hashmap1_init2(gs, num_gs);

    gs_base = n_mulmod_shoup(a, bs[num_bs - 1], bs_shoup[0], n);
    gs_shoup = n_mulmod_precomp_shoup(gs_pow, n);
    gs_pow = gs_base;

    for (i = 1; i < num_gs; i++)
    {
        hashmap1_insert(gs_pow, (void *) (i*num_bs), gs);
        gs_pow = n_mulmod_shoup(gs_base, gs_pow, gs_shoup, n);
    }
    hashmap1_insert(gs_pow, (void *) (i*num_bs), gs);

    for (i = 0; i < num; i++)
    {

        if (hashmap1_find((void **) &val, b[i], gs))
        {
            r[i] = val;

            continue;
        }

        for (j = 0; j < num_bs; j++)
        {
            key = n_mulmod_shoup(bs[j], b[i], bs_shoup[j], n);
            if (hashmap1_find((void **) &val, key, gs))
            {
                r[i] = val - j - 1;

                break;
            }
        }
    }

    hashmap1_clear(gs);
}

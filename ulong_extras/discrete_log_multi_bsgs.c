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
    ulong num_gs, num_bs, gs_base, gs_pow, gs_shoup, gs_gap;
    ulong *bs, *bs_shoup;
    ulong val, key;
    hashmap1_t gs;

    if (num == 0)
        return;

    if (ord == 0)
        ord = n - 1;

    /* number of giant steps */
    num_gs = ceil(sqrt((double)ord * num));

    /* number of baby steps = ceil(ord/num_gs) - 1*/
    num_bs = (ord - 1) / num_gs;

    /* giant steps gap = the power of a corresponding to first giant step */
    gs_gap = num_bs + 1;

    hashmap1_init2(gs, num_gs);

    /* fill in baby steps table */
    if (num_bs)
    {
        bs = flint_malloc(num_bs * sizeof(ulong));
        bs_shoup = flint_malloc(num_bs * sizeof(ulong));
        bs[0] = a;
        bs_shoup[0] = n_mulmod_precomp_shoup(a, n);

        for (i = 1; i < num_bs; i++)
        {
            bs[i] = n_mulmod_shoup(a, bs[i - 1], bs_shoup[0], n);
            bs_shoup[i] = n_mulmod_precomp_shoup(bs[i], n);
        }

        gs_base = n_mulmod_shoup(a, bs[num_bs - 1], bs_shoup[0], n);
        gs_shoup = n_mulmod_precomp_shoup(gs_base, n);
    }
    else
    {
        bs = NULL;
        bs_shoup = NULL;
        gs_base = a;
        gs_shoup = n_mulmod_precomp_shoup(a, n);
    }

    /* fill in giant steps table */
    gs_pow = gs_base;
    for (i = 1; i < num_gs; i++)
    {
        hashmap1_insert(gs_pow, (void *) (i * gs_gap), gs);
        gs_pow = n_mulmod_shoup(gs_base, gs_pow, gs_shoup, n);
    }
    hashmap1_insert(gs_pow, (void *) (i * gs_gap), gs);

    /* compute baby steps for each element and perform lookups */
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

        FLINT_ASSERT (j < num_bs);
    }

    if (num_bs)
    {
        flint_free(bs);
        flint_free(bs_shoup);
    }

    hashmap1_clear(gs);
}

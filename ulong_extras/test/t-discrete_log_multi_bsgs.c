/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"


int main(void)
{
    int i, j;
    ulong *pows, *logs;
    const ulong maxnum = 100;
    FLINT_TEST_INIT(state);

    flint_printf("discrete_log_multi_bsgs....");
    fflush(stdout);

    flint_randinit(state);

    pows = flint_malloc(maxnum * sizeof *pows);
    logs = flint_malloc(maxnum * sizeof *pows);

    /* tests with big primes and unknown order */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong p, base, num, check = 0;
        double pinv;

        p = n_randprime(state, n_randint(state, 16) + 10, 1);
        pinv = n_precompute_inverse(p);
        num = n_randint(state, maxnum + 1);

        base = n_urandint(state, p);

        for (j = 0; j < num; ++j)
        {
            pows[j] = n_powmod_precomp(base, n_randint(state, p - 1), p, pinv);
        }

        n_discrete_log_multi_bsgs(logs, pows, base, num, 0, p);

        for (j = 0; j < num; ++j)
        {
            check = n_powmod_precomp(base, logs[j], p, pinv);
            if (check != pows[j])
                break;
        }

        if (j < num)
        {
            flint_printf("FAIL1:\n");
            flint_printf("array index %d of %wu\n", j, num);
            flint_printf("%wu ** %wu mod %wu == %wu != %wu\n",
                    base, logs[j], p, check, pows[j]);
            abort();
        }
    }

    /* tests with small primes and known order */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong p, base, num, check = 0;
        ulong ord, bpow;
        double pinv;

        p = n_randprime(state, 10, 1);
        pinv = n_precompute_inverse(p);
        num = n_randint(state, maxnum + 1);

        base = n_urandint(state, p - 1) + 1;
        bpow = base;
        for (ord = 1; bpow != 1; ++ord)
        {
            bpow = n_mulmod_precomp(bpow, base, p, pinv);
        }

        for (j = 0; j < num; ++j)
        {
            pows[j] = n_powmod_precomp(base, n_randint(state, ord), p, pinv);
        }

        n_discrete_log_multi_bsgs(logs, pows, base, num, ord, p);

        for (j = 0; j < num; ++j)
        {
            check = n_powmod_precomp(base, logs[j], p, pinv);
            if (check != pows[j])
                break;
        }

        if (j < num)
        {
            flint_printf("FAIL2:\n");
            flint_printf("array index %d of %wu\n", j, num);
            flint_printf("%wu ** %wu mod %wu == %wu != %wu\n",
                    base, logs[j], p, check, pows[j]);
            abort();
        }
    }

    flint_free(pows);
    flint_free(logs);
    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

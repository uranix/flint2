/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

#include "profiler.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    

    flint_printf("get_cfrac....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_t x, r;
        fmpz *c1, *c2;
        slong n1, n2, bound;

        fmpq_init(x);
        fmpq_init(r);

        fmpq_randtest(x, state, 1 + n_randint(state, 2000));
        bound = fmpq_cfrac_bound(x);

        c1 = _fmpz_vec_init(bound);
        c2 = _fmpz_vec_init(bound);

        n1 = fmpq_get_cfrac(c1, r, x, bound);

        if (!fmpq_is_zero(r))
        {
            flint_printf("FAIL: expected zero remainder\n");
            abort();
        }

        /* Test chaining */
        n2 = 0;
        while (1)
        {
            n2 += fmpq_get_cfrac(c2 + n2, x, x, n_randint(state, 10));
            if (fmpq_is_zero(x))
                break;
            fmpq_inv(x, x);
        }

        if (n1 != n2)
        {
            flint_printf("FAIL: i = %wd, n1 = %wd, n2 = %wd\n", i, n1, n2);
            abort();
        }

        if (!_fmpz_vec_equal(c1, c2, n1))
        {
            flint_printf("FAIL: vectors not equal\n");
            _fmpz_vec_print(c1, n1); flint_printf("\n");
            _fmpz_vec_print(c2, n2); flint_printf("\n");
            abort();
        }

        _fmpz_vec_clear(c1, bound);
        _fmpz_vec_clear(c2, bound);
        fmpq_clear(x);
        fmpq_clear(r);
    }


    for (i = 0; i <= 20; i++)
    {
        fmpq_t x, r;
        fmpz *c1, *c2;
        slong n1, n2, bound;
        slong gcd_time;
timeit_t timer;

        fmpq_init(x);
        fmpq_init(r);

        fmpz_fib_ui(fmpq_numref(x), (1 << i) + 1);
        fmpz_fib_ui(fmpq_denref(x), (1 << i));

flint_printf("--- fib i = %wd (numerator bits = %wu) ---\n", i, fmpz_bits(fmpq_numref(x)));

timeit_start(timer);
        fmpq_canonicalise(x);
timeit_stop(timer);
gcd_time = timer->wall;
flint_printf("gcd: %wd\n", timer->wall);

        bound = fmpq_cfrac_bound(x);
        c1 = _fmpz_vec_init(bound);
        c2 = _fmpz_vec_init(bound);

timeit_start(timer);
        n1 = fmpq_get_cfrac(c1, r, x, bound);
timeit_stop(timer);
        if (!fmpq_is_zero(r))
        {
            flint_printf("FAIL1: expected zero remainder\n");
            abort();
        }
gcd_time = FLINT_MAX(1, gcd_time);
flint_printf("new: %wd  (new/gcd: %f)\n", timer->wall, (double)(timer->wall)/(double)(gcd_time));


timeit_start(timer);
        n2 = fmpq_get_cfracOLD(c2, r, x, bound);
timeit_stop(timer);
        if (!fmpq_is_zero(r))
        {
            flint_printf("FAIL2: expected zero remainder\n");
            abort();
        }

flint_printf("old: %wd\n", timer->wall);

        if (n1 != n2)
        {
            flint_printf("FAIL3: n1 != n2\n");
            abort();
        }

        _fmpz_vec_clear(c1, bound);
        _fmpz_vec_clear(c2, bound);
        fmpq_clear(x);
        fmpq_clear(r);

    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

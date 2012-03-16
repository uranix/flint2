/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_dmpoly.h"

int main(void)
{
    flint_rand_t state;
    long iter;

    printf("set_coeff_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        nmod_dmpoly_t A;
        mp_limb_t mod, a, b;
        long index[10];
        long i, j;
        int vars;

        mod = n_randtest_prime(state, 0);
        vars = 1 + n_randint(state, 4);

        nmod_dmpoly_init(A, vars, mod);
        nmod_dmpoly_randtest(A, state, n_randint(state, 10));

        for (i = 0; i < 10; i++)
        {
            for (j = 0; j < vars; j++)
                index[j] = n_randint(state, 10);

            a = n_randtest(state) % mod;
            nmod_dmpoly_set_coeff_ui(A, index, a);
            b = nmod_dmpoly_get_coeff_ui(A, index);

            if (a != b)
            {
                printf("FAIL\n");
                printf("A:\n"); nmod_dmpoly_print(A); printf("\n\n");
                printf("mod: %lu\n", mod);
                printf("a = %lu, b = %lu\n", a, b);
                abort();
            }
        }

        nmod_dmpoly_clear(A);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}

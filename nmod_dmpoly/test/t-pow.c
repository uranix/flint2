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

    printf("pow....");
    fflush(stdout);

    flint_randinit(state);

    /* Check evaluation homomorphism */
    for (iter = 0; iter < 10000; iter++)
    {
        nmod_dmpoly_t A, B;
        mp_limb_t mod, va, vb;
        mp_limb_t points[10];
        ulong exp;
        long len, i;
        int vars;

        mod = n_randtest_prime(state, 0);
        len = n_randint(state, 5);
        vars = 1 + n_randint(state, 3);
        exp = n_randint(state, 10);

        for (i = 0; i < vars; i++)
            points[i] = n_randlimb(state) % mod;

        nmod_dmpoly_init(A, vars, mod);
        nmod_dmpoly_init(B, vars, mod);

        nmod_dmpoly_randtest(A, state, len);
        if (n_randint(state, 2))
            nmod_dmpoly_randtest(B, state, len);

        nmod_dmpoly_pow(B, A, exp);

        va = nmod_dmpoly_evaluate_nmod(A, points);
        vb = nmod_dmpoly_evaluate_nmod(B, points);

        if (vb != n_powmod2_preinv(va, exp, mod, n_preinvert_limb(mod)))
        {
            printf("FAIL\n");
            printf("A:\n"); nmod_dmpoly_print(A); printf("\n\n");
            printf("B:\n"); nmod_dmpoly_print(B); printf("\n\n");
            printf("exp: %lu\n", exp);
            printf("mod: %lu\n", mod);
            printf("va = %lu, vb = %lu\n", va, vb);
            abort();
        }

        nmod_dmpoly_clear(A);
        nmod_dmpoly_clear(B);
    }

    /* Check aliasing */
    for (iter = 0; iter < 1000; iter++)
    {
        nmod_dmpoly_t A, B;
        mp_limb_t mod;
        ulong exp;
        long len;
        int vars;

        mod = n_randtest_prime(state, 0);
        len = n_randint(state, 5);
        vars = 1 + n_randint(state, 3);
        exp = n_randint(state, 10);

        nmod_dmpoly_init(A, vars, mod);
        nmod_dmpoly_init(B, vars, mod);

        nmod_dmpoly_randtest(A, state, len);
        if (n_randint(state, 2))
            nmod_dmpoly_randtest(B, state, len);

        nmod_dmpoly_pow(B, A, exp);
        nmod_dmpoly_pow(A, A, exp);

        if (!nmod_dmpoly_equal(A, B))
        {
            printf("FAIL\n");
            printf("A:\n"); nmod_dmpoly_print(A); printf("\n\n");
            printf("B:\n"); nmod_dmpoly_print(B); printf("\n\n");
            printf("exp: %lu\n", exp);
            printf("mod: %lu\n", mod);
            abort();
        }

        nmod_dmpoly_clear(A);
        nmod_dmpoly_clear(B);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}

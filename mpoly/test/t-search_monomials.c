/*
    Copyright (C) 2016 William Hart

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
#include "mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    ulong a[]={289,278,249,224,217,215,206,192,165,136,106,89,85,83,61,58,34,21,0};
    ulong b[]={111,87,83,56,53,41,33,24,0};
    ulong e[]={0,0,0,0}
    slong escore;

    flint_printf("search_monomials....");
    fflush(stdout);

    fmpz_mpoly_search_monomials(e,&escore, 34, 34, a, 19, b, 9, FLINT_BITS, 1, 0, 0);

    flint_printf("e: %wu",e[0]);
    flint_printf("escore: %wd",escore);



    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}


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
    ulong b[4]={8,8,0,0};
    ulong a[4]={5,0,1,4};
    ulong x[4]={59,59,59,59};
    ulong m[4]={96,97,98,99};


    flint_printf("average....");
    fflush(stdout);

    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
    }

    _mpoly_average_monomial_degrevlex(x, a, b, m, 4);

   flint_printf("{");
   for (i=0; i<4; i++)
   {
      flint_printf("%wd",x[i]);
      if(i+1<4)
         flint_printf(",");      
   }
   flint_printf("}\n");


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}


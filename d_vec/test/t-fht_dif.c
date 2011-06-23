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

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "d_vec.h"

int main(void)
{
   int result;
   long i, j;
   flint_rand_t state;
   
   printf("fht_dif_rec/iter....");
   fflush(stdout);

   flint_randinit(state);
   
   for (i = 0; i < 24; i++)
   {
      long n = (1L<<i);

      double * stab = _d_vec_init((5*n)/16);
      double * vec = _d_vec_init(n);
      double * vec2 = _d_vec_init(n);
      _d_vec_randbits(vec, state, n, 53);
      for (j = 0; j < n; j++)
          vec2[j] = vec[j];

      _d_vec_compute_stab(stab, n);
      _d_vec_fht_dif_rec(vec, stab, 1, n);
      _d_vec_fht_dif_iter(vec2, stab, 1, n);

      result = 1;
      for (j = 0; j < n; j++)
      {
          if (vec[j] != vec2[j])
          {
              result = 0;
              break;
          }
      }

      if (!result)
      {
         printf("FAIL:\n");
         printf("n = %ld\n", n);
         printf("failure at j = %ld\n", j);
         printf("vec[j] = %f, vec2[j] = %f\n", vec[j], vec2[j]);
         abort();
      }

      _d_vec_free(vec);
      _d_vec_free(vec2);
      _d_vec_free(stab);
   }

   flint_randclear(state);

   printf("PASS\n");
   return 0;
}

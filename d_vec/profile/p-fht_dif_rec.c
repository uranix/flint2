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

    Copyright 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "d_vec.h"

typedef struct
{
   long n;
} info_t;

void sample1(void * arg, ulong count)
{
   info_t * info = (info_t *) arg;
   long n = info->n, i, j;
   long scale;

   flint_rand_t state;
   flint_randinit(state);
     
   double * stab = _d_vec_init((5*n)/16);

   scale = 1;
   if (n < 100000) scale = 10;
   if (n < 10000) scale = 100;
   if (n < 100) scale = 1000;
      
   for (i = 0; i < count; i++)
   {
	  prof_start();
      for (j = 0; j < scale; j++)
          _d_vec_compute_stab(stab, n);
	  prof_stop();
      if (stab[n/4 - 1] == 0.1f) abort();
   }
  
   _d_vec_free(stab);
   flint_randclear(state);
}

void sample2(void * arg, ulong count)
{
   info_t * info = (info_t *) arg;
   long n = info->n, i, j;
   long scale;

   flint_rand_t state;
   flint_randinit(state);
     
   double * stab = _d_vec_init((5*n)/16);
   double * vec = _d_vec_init(n);
   _d_vec_randbits(vec, state, n, 53);

   scale = 1;
   if (n < 100000) scale = 10;
   if (n < 10000) scale = 100;
   if (n < 100) scale = 1000;
      
   for (i = 0; i < count; i++)
   {
	  prof_start();
      for (j = 0; j < scale; j++)
      {
          _d_vec_compute_stab(stab, n);
          _d_vec_fht_dif_rec(vec, stab, 1, n);
      }
	  prof_stop();
      if (stab[n/4 - 1] == 0.1f) abort();
   }
  
   _d_vec_free(vec);
   _d_vec_free(stab);
   flint_randclear(state);
}

int main(void)
{
   double min, max, min2, max2;
   info_t info;
   long k, scale;

   for (k = 2; k <= 30; k++)
   {
      info.n = 1<<k;

      scale = 1;
      if (info.n < 100000) scale = 10;
      if (info.n < 10000) scale = 100;
      if (info.n < 100) scale = 1000;
   
      prof_repeat(&min, &max, sample1, (void *) &info);
         
      printf("length %d, min %.3g ms, max %.3g ms\n", 
           info.n,
		   ((min/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((max/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0
	     );
     
      prof_repeat(&min2, &max2, sample2, (void *) &info);
         
      printf("       min %.3g ms, max %.3g ms, norm %.3g\n", 
           (((min2 - min)/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           (((max2 - max)/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0,
           ((((min2 - min)/(double)FLINT_CLOCK_SCALE_FACTOR)/scale)/2400000.0)
              *500000.0/info.n/FLINT_BIT_COUNT(info.n)
	     );
   }

   return 0;
}

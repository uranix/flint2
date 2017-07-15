/*
    Copyright (C) 2017

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "mpoly.h"


void _mpoly_average_monomial_lex(ulong * x, ulong * A, ulong * B,
                                                          ulong * max, slong n)
{
   int pref = 0;
   slong i, j;
   ulong * a, * b, * aedge, * bedge;
   ulong adiff, bdiff;
   TMP_INIT;

   TMP_START;
   a = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   b = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   aedge = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   bedge = (ulong *) TMP_ALLOC(n*sizeof(ulong));

   for (i = 0; i < n; i++)
   {
      a[i] = A[i];
      b[i] = B[i];
   }

   while (1)
   {

flint_printf("\nb: ");
for (i = 0; i < n; i++)
flint_printf(" %wu",b[i]);

flint_printf("\na: ");
for (i = 0; i < n; i++)
flint_printf(" %wu",a[i]);

flint_printf("\n");


      for (i = 0; i < n - 1; i++)
      {
         aedge[i] = a[i];
         bedge[i] = b[i];
         x[i] = a[i];

         FLINT_ASSERT(b[i] >= a[i]);

         if (b[i] - a[i] > 1)
         {
            x[i] = (a[i] + b[i] + 1)/2;
            for (j = i + 1; j < n; j++)
               x[j] = 0;
            TMP_END;
            return;

         } else if (b[i] - a[i] == 1)
         {
            adiff = bdiff = 0;
            for (j = i + 1; j < n; j++)
            {
               adiff |= a[j] - max[j];
               bdiff |= b[j];
               aedge[j] = max[j];
               bedge[j] = 0;
            }

            if (adiff == 0 && bdiff == 0)
            {
               for (j = 0; j < n; j++)
                  x[j] = (pref == 0 ? aedge : bedge)[j];
               TMP_END;
               return;

            } else if (adiff == 0 && bdiff != 0)
            {
               for (j = 0; j < n; j++)
                  a[j] = bedge[j];
               pref = 0;
               break;

            } else if (adiff != 0 && bdiff == 0)
            {
               for (j = 0; j < n; j++)
                  b[j] = aedge[j];
               pref = 1;
               break;

            } else
            {
               for (j = 0; j < n; j++)
                  x[j] = aedge[j];
               TMP_END;
               return;
            }
         }
      }
   }
}




void _mpoly_average_monomial_deglex(ulong * x, ulong * A, ulong * B,
                                                          ulong * max, slong n)
{
   int pref = 0;
   slong i, j;
   ulong * a, * b, * aedge, * bedge;
   ulong adiff, bdiff, atot, btot;
   TMP_INIT;

   TMP_START;
   a = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   b = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   aedge = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   bedge = (ulong *) TMP_ALLOC(n*sizeof(ulong));

   for (i = 0; i < n; i++)
   {
      a[i] = A[i];
      b[i] = B[i];
   }

   while (1)
   {
/*
flint_printf("\nb: ");
for (i = 0; i < n; i++)
flint_printf(" %wu",b[i]);

flint_printf("\na: ");
for (i = 0; i < n; i++)
flint_printf(" %wu",a[i]);

flint_printf("\n");
*/

      atot = btot = 0;
      for (i = 0; i < n - 1; i++)
      {
         aedge[i] = a[i];
         bedge[i] = b[i];
         x[i] = a[i];

         if (i != 0)
         {
            atot += a[i];
            btot += b[i];
         }

         FLINT_ASSERT(b[i] >= a[i]);

         if (b[i] - a[i] > 1)
         {
            x[i] = (a[i] + b[i] + 1)/2;
            if (i != 0)
               atot += x[i] - aedge[i];
            for (j = i + 1; j < n - 1; j++)
               x[j] = 0;
            x[n - 1] = x[0] - atot;
            TMP_END;
            return;

         } else if (b[i] - a[i] == 1)
         {
            adiff = bdiff = 0;
            for (j = i + 1; j < n - 1; j++)
            {
               adiff |= a[j + 1];
               bdiff |= b[j];
               aedge[j + 1] = 0;
               bedge[j] = 0;
            }
            bedge[n - 1] = bedge[0] - btot;
            aedge[i + 1] = aedge[0] - atot;

            if (adiff == 0 && bdiff == 0)
            {
               for (j = 0; j < n; j++)
                  x[j] = (pref == 0 ? aedge : bedge)[j];
               TMP_END;
               return;

            } else if (adiff == 0 && bdiff != 0)
            {
               for (j = 0; j < n; j++)
                  a[j] = bedge[j];
               pref = 0;
               break;

            } else if (adiff != 0 && bdiff == 0)
            {
               for (j = 0; j < n; j++)
                  b[j] = aedge[j];
               pref = 1;
               break;

            } else
            {
               for (j = 0; j < n; j++)
                  x[j] = aedge[j];
               TMP_END;
               return;
            }
         }
      }
   }
}



void _mpoly_average_monomial_degrevlex(ulong * x, ulong * A, ulong * B,
                                                          ulong * max, slong n)
{
   int pref = 0;
   slong i, j;
   ulong * a, * b, * aedge, * bedge;
   ulong cmp, adiff, bdiff, atot, btot;
   TMP_INIT;

   TMP_START;
   a = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   b = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   aedge = (ulong *) TMP_ALLOC(n*sizeof(ulong));
   bedge = (ulong *) TMP_ALLOC(n*sizeof(ulong));

   for (i = 0; i < n; i++)
   {
      a[i] = A[i];
      b[i] = B[i];
   }

   while (1)
   {
/*
flint_printf("\nb: ");
for (i = 0; i < n; i++)
flint_printf(" %wu",b[i]);

flint_printf("\na: ");
for (i = 0; i < n; i++)
flint_printf(" %wu",a[i]);

flint_printf("\n");
*/

      atot = btot = 0;
      for (i = 0; i < n - 1; i++)
      {
         aedge[i] = a[i];
         bedge[i] = b[i];
         x[i] = a[i];

         if (i != 0)
         {
            atot += a[i];
            btot += b[i];
         }

         FLINT_ASSERT((i==0 && b[i] >= a[i])||(i!=0 && b[i] <= a[i]));

         cmp = b[i] - a[i];
         if (i != 0)
            cmp = -cmp;

         if (cmp > 1)
         {
            x[i] = (a[i] + b[i])/2;
            if (i != 0)
               atot += x[i] - aedge[i];
            for (j = i + 1; j < n - 1; j++)
               x[j] = 0;
            x[n - 1] = x[0] - atot;
            TMP_END;
            return;

         } else if (cmp == 1)
         {
            adiff = bdiff = 0;
            for (j = i + 1; j < n - 1; j++)
            {
               adiff |= a[j];
               bdiff |= b[j + 1];
               aedge[j] = 0;
               bedge[j + 1] = 0;
            }
            bedge[i + 1] = bedge[0] - btot;
            aedge[n - 1] = aedge[0] - atot;

            if (adiff == 0 && bdiff == 0)
            {
flint_printf("here1\n");
               for (j = 0; j < n; j++)
                  x[j] = (pref == 0 ? aedge : bedge)[j];
               TMP_END;
               return;

            } else if (adiff == 0 && bdiff != 0)
            {
flint_printf("here2\n");
               for (j = 0; j < n; j++)
                  a[j] = bedge[j];
               pref = 0;
               break;

            } else if (adiff != 0 && bdiff == 0)
            {
flint_printf("here3\n");
               for (j = 0; j < n; j++)
                  b[j] = aedge[j];
               pref = 1;
               break;

            } else
            {
flint_printf("here4\n");
               for (j = 0; j < n; j++)
                  x[j] = aedge[j];
               TMP_END;
               return;
            }
         }
      }
   }
}


/*
   average the monomials a and b into the monomial x
   assume that a <= b
   assume that the fields of a and b are both <= fields of the unpacked abmax
   x is guaranteed to satisfy a<=x<=b if a<=b
   return is nonzero iff a != x
*/
ulong mpoly_average_monomial(ulong * x, ulong * a, ulong * b,
              ulong * abmax, slong bits, slong nfields, int deg, int rev)
{
   slong i;
   ulong res = 0;
   ulong * unp_x, * unp_a, * unp_b;

   unp_x = flint_malloc(nfields*sizeof(ulong));
   unp_a = flint_malloc(nfields*sizeof(ulong));
   unp_b = flint_malloc(nfields*sizeof(ulong));

/*
   unpack both exponents
*/

   mpoly_unpack_monomials_noalloc(unp_a, FLINT_BITS, a, 1, nfields, bits);
   mpoly_unpack_monomials_noalloc(unp_b, FLINT_BITS, b, 1, nfields, bits);

   if (deg==0 && rev==0)
      _mpoly_average_monomial_lex(unp_x, unp_a, unp_b, abmax, nfields);
   else if (deg!=0 && rev==0)
      _mpoly_average_monomial_deglex(unp_x, unp_a, unp_b, abmax, nfields);
   else if (deg!=0 && rev!=0)
      _mpoly_average_monomial_degrevlex(unp_x, unp_a, unp_b, abmax, nfields);

   for (i = 0; i < nfields; i++)
      res |= unp_x[i] ^ unp_a[i];
/*
   hacky - we just want to pack all fields from FLINT_BITS to bits
*/
   mpoly_get_monomial(a, unp_a, bits, nfields, 0, 0);

   flint_free(unp_b);
   flint_free(unp_b);
   flint_free(unp_x);

   return res;
}


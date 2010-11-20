#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "generic.h"

void test1(void)
{
   long i;
   
   vec_t r, a;
   obj_t c;

   r.typ = DOUBLE;
   a.typ = DOUBLE;
   r.len = 1000;
   a.len = 1000;
   r.arr = malloc(1000*sizeof(double));
   a.arr = malloc(1000*sizeof(double));

   c.typ = DOUBLE;
   c.val.dbl = 1.000001;

   for (i = 0; i < 1000000; i++)
      vec_scalar_mul(&r, &a, &c);
}

void test2(void)
{
   long i;
   
   double * r, * a;
   double c;

   r = malloc(1000*sizeof(double));
   a = malloc(1000*sizeof(double));

   c = 1.000001;

   for (i = 0; i < 1000000; i++)
      vecmul_d(r, a, 1000, c);
}

int main(void)
{
   test1();

   return 0;
}

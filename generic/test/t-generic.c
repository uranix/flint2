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

void test3(void)
{
   long i;
   
   fmpz_poly_t * r, * a;
   fmpz_poly_t c;
   flint_rand_t state;

   flint_randinit(state);

   r = malloc(1000*sizeof(fmpz_poly_t));
   a = malloc(1000*sizeof(fmpz_poly_t));

   for (i = 0; i < 1000; i++)
   {
      fmpz_poly_init(r[i]);
	  fmpz_poly_init(a[i]);
	  fmpz_poly_randtest(a[i], state, 10, 50);     
   }

   fmpz_poly_init(c);
   fmpz_poly_randtest(c, state, 10, 50);     

   for (i = 0; i < 1000; i++)
      vecmul_fmpz_poly(r, a, 1000, c);

   flint_randclear(state);

}

void test4(void)
{
   long i;
   
   flint_rand_t state;

   vec_t r, a;
   obj_t c;

   flint_randinit(state);

   r.typ = FMPZ_POLY;
   a.typ = FMPZ_POLY;
   r.len = 1000;
   a.len = 1000;
   r.arr = malloc(1000*sizeof(fmpz_poly_t));
   a.arr = malloc(1000*sizeof(fmpz_poly_t));

   for (i = 0; i < 1000; i++)
   {
      fmpz_poly_init(((fmpz_poly_t *)r.arr)[i]);
	  fmpz_poly_init(((fmpz_poly_t *)a.arr)[i]);
	  fmpz_poly_randtest(((fmpz_poly_t *)a.arr)[i], state, 10, 50);     
   }

   c.typ = FMPZ_POLY;
   fmpz_poly_init(c.val.f_poly);
   fmpz_poly_randtest(c.val.f_poly, state, 10, 50);     
  
   for (i = 0; i < 1000; i++)
	  vec_scalar_mul(&r, &a, &c);

   flint_randclear(state);
}

int main(void)
{
   test3();

   return 0;
}

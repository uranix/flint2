#include <mpir.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "generic.h"

void vecmul_d(double * r, double * a, ulong len, double c)
{
   long i;
   for (i = 0; i < len; i++)
      r[i] = a[i] * c;
}

void vecmul_fmpz_poly(fmpz_poly_t * r, fmpz_poly_t * a, ulong len, fmpz_poly_t c)
{
   long i;
   for (i = 0; i < len; i++)
      fmpz_poly_mul(r[i], a[i], c);
}

void vecmul_ui(mp_limb_t * r, mp_limb_t * a, ulong len, mp_limb_t c)
{
   long i;
   for (i = 0; i < len; i++)
      r[i] = a[i] * c;
}

void mul_d(void * r, void * a, void * c)
{
	*(double *)r = *(double *)a * *(double *)c;
}

void mul_ui(void * r, void * a, void * c)
{
	*(mp_limb_t *)r = *(mp_limb_t *)a * *(mp_limb_t *)c;
}

void mul_fmpz_poly(void * r, void * a, void * c)
{
	fmpz_poly_mul(*(fmpz_poly_t *)r, *(fmpz_poly_t *)a, *(fmpz_poly_t *)c);
}

typedef void (* genmul_t)(void *, void *, void *);

void vec_scalar_mul(vec_t * r, vec_t * a, obj_t * c)
{
   ulong len = a->len;
   long i;
   size_t size = 0;
   genmul_t fn;

   switch (c->typ)
   {
   case DOUBLE:
      size = sizeof(double);
	  fn = mul_d;
	  break;
   case ULONG:
      size = sizeof(ulong);
	  fn = mul_ui;
	  break;
   case FMPZ_POLY:
      size = sizeof(fmpz_poly_t);
	  fn = mul_fmpz_poly;
	  break;
   default:
      fn = NULL;
   }

   for (i = 0; i < len; i++)
      fn(r->arr + i*size, a->arr + i*size, &(c->val));
}
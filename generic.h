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

    Copyright (C) 2010 William Hart
 
******************************************************************************/

#ifndef GENERIC_H
#define GENERIC_H

#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz_poly.h"

typedef enum type_t
{
   DOUBLE, ULONG, FMPZ_POLY
} type_t;

typedef struct obj_t
{
   type_t typ;
   
   union
   {
      double dbl;
	  ulong ui;
	  fmpz_poly_t f_poly;
   } val;
} obj_t;

typedef struct vec_t
{
   type_t typ;
   ulong len;

   char * arr;
} vec_t;

void vecmul_d(double * r, double * a, ulong len, double c);

void vecmul_fmpz_poly(fmpz_poly_t * r, fmpz_poly_t * a, ulong len, fmpz_poly_t c);

void vec_scalar_mul(vec_t * r, vec_t * a, obj_t * c);

#endif


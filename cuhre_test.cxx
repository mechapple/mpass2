/*
 * cuhre_test.cxx
 * 
 * Copyright 2022 Anirban <anirban@anirban-ThinkPad-T430>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <iostream>
#include <cmath>
#include "cuba.h"

#define NDIM 2
#define NCOMP 12
#define NVEC 1
#define EPSREL 1e-6
#define EPSABS 0
#define VERBOSE 0
#define LAST 4

#define MINEVAL 0
#define MAXEVAL 50000
#define STATEFILE NULL
#define SPIN NULL
#define KEY 0

static int f12(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  for(int i=0;i<12;i++)
    ff[i] = 1.0/(1.0+i*xx[0]*xx[0]);
  
  return 0;
}

int main(int argc, char **argv)
{
  int comp, nregions, neval, fail;
  double integral[12], error[12], prob[12];
    
  Cuhre(NDIM, NCOMP, f12, NULL, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
    
  for(int i=0;i<12;i++) printf("%lf ",integral[i]); printf("\n");

	return 0;
}


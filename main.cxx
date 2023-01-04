/*
 * main.cxx
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
//g++ -O3 -g -std=c++11 -frounding-math main.cxx -o main -lm -lgsl -lgslcblas -lcuba -lnlopt

#include "headers.h"

int main(int argc, char **argv)
{
  fibnetwork fn;
  fn.readfile(argv[1]);
  //fn.printdata();
  fn.printlammps_cps((char*) "Output_cps.lammpstrj",(char*) "w");
  fn.printlammps((char*) "Output.lammpstrj",(char*) "w", 20);

  //clock_t begin1,end1;
  //begin1 = clock();
  //fn.compute_ef();
  
  //end1 = clock();
  //double time1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
  //printf ("time = %lf\n",time1);
  
  
  if(1) loop(i,20) {
	fn.minimize();
  //fn.integrate_runge_kutta_4(10000,1e-4);
	integrate(fn,1e2,2e-3);
  }
  //fn.printlammps_cps((char*) "Output_cps.lammpstrj",(char*) "a");

  return 0;
}

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
// g++ -O3 -g -frounding-math main.cxx -o main -lm -lgsl -lgslcblas -lcuba -lnlopt

#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <chrono>

#include <nlopt.h>
#include <dlib/optimization.h>


// INCLUDE GSL HEADERS
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#define De 0.015
#define re 0.3
#define am 10.0 //2*alpha
#define bm 5.0 //alpha

#define E 3000
#define R 0.3
#define gamma 0.015
#define L0 21
#define u0 0.3
#define x1 20
#define u1 5

#define loop(x,n) for(int x = 0; x < n; ++x)
#define loop2(x,a,b) for(int x = a; x < b; ++x)

int NPOINTS = 8;

struct Vector3
{
  double comp[3];

  Vector3(double x, double y, double z)
  {
    comp[0] = x;
    comp[1] = y;
    comp[2] = z;
  }

  Vector3()
  {
    comp[0] = 0;
    comp[1] = 0;
    comp[2] = 0;
  }
};

double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
  return (x[0] - x[1] - x1 * 0.1);
}

double VdV(unsigned n, const double *x, double *grad, void *fn_data)
{
  double l1 = x[0], l2 = x[1], u2 = x[2];
  double EA_2L = E * M_PI * R * R / (2 * L0);
  double EI = 0.25 * E * M_PI * R * R * R * R;

  double a = 2 * (u2 - u0) / (l2 - l1) - (u1 - u2) / (x1 - l2);
  double b = ((u1 - u2) / (x1 - l2) - (u2 - u0) / (l2 - l1)) / (l2 - l1);

  double L13 = l1 + sqrt((u1 - u2) * (u1 - u2) + (x1 - l2) * (x1 - l2));
  double t2 = 2 * b * (l2 - l1);
  double L2 = ((a + t2) * sqrt(1 + (a + t2) * (a + t2)) + asinh(a + t2) - a * sqrt(1 + a * a) - asinh(a)) / (4 * b);
  double Lt = L13 + L2;

  double Ua = EA_2L * (Lt - L0) * (Lt - L0);

  double B2 = (a + t2) * (2 * (a + t2) * (a + t2) + 3) / (3 * pow(1 + (a + t2) * (a + t2), 1.5)) - a * (2 * a * a + 3) / (3 * pow(1 + a * a, 1.5));
  double Ub = EI * 2 * b * b * B2;
  double Uc = -gamma * l1;

  return (Ua + Ub + Uc);
}

#define _DMAX 1000000

typedef dlib::matrix<double,0,1> column_vector;

struct Point
{
    double x;
    double y;
};

static Point points1[] =
{
    { 254 , 102 },
    { 226 , 63  },
    { 185 , 49  },
    { 146 , 74  },
    { 142 , 119 },
    { 117 , 169 },
    { 86  , 214 },
    { 40  , 200 },
};

//
// cubicBezier
//
// p1 - start point
// c1 - first control point
// c2 - second control point
// p2 - end point
//
double cubicBezier(double p1, double c1, double c2, double p2, double t)
{
    double s = (1 - t);

    double v = 0;
    v = v + (1 * p1 * s * s * s);
    v = v + (3 * c1 * s * s * t);
    v = v + (3 * c2 * s * t * t);
    v = v + (1 * p2 * t * t * t);

    return v;        
};

Point cubicBezier(double p1x, double p1y, double c1x, double c1y,
                  double c2x, double c2y, double p2x, double p2y,
                  double t)
{
    Point pt;
    pt.x = cubicBezier(p1x, c1x, c2x, p2x, t);
    pt.y = cubicBezier(p1y, c1y, c2y, p2y, t);
    return pt;
}

// Any distance function can be used for optimisation.  This one, where we want
// to find the least-squares is most common. 
double dist(double dx1, double dy1, double dx2, double dy2)
{
    double x = dx2 - dx1;
    double y = dy2 - dy1;

    return (x * x) + (y * y);
}

// This is the function that the optimiser calls repeatedly with different
// parameters, attempting to get the lowest return value it can.
double rateCurve(const column_vector& params)
{
    double p1x = points1[0].x;
    double p1y = points1[0].y;
    double c1x = params(0,0);
    double c1y = params(1,0);
    double c2x = params(2,0);
    double c2y = params(3,0);
    double p2x = points1[NPOINTS-1].x;
    double p2y = points1[NPOINTS-1].y;

    double distances = 0;

    for (Point target : points1)
    {
        double distance = _DMAX;

        for (double t = 0; t <= 1; t += 0.02)
        {
            Point pt = cubicBezier( p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y, t);

            double dCandidate = dist(pt.x, pt.y, target.x, target.y);

            distance = std::min(distance, dCandidate);
        }

        distances += distance;
    }

    // Thats the curve-fitting done.  Now incentivise slightly smoother curves.
    double p1c1 = dist(p1x, p1y, c1x, c1y);
    double p2c2 = dist(p2x, p2y, c2x, c2y);

    return distances + pow(p1c1, 0.6) + pow(p2c2, 0.6);
}

int main(int argc, char **argv)
{
  double *x = new double[3];
  double *df = new double[3];
  // loop(i,np) loop(j,3) x[i*3+j] = CPx[i].comp[j];

  nlopt_opt opt;

  // establish sizes
  unsigned n = 3;    // number of decision variables
  unsigned m_in = 0; // number of inequality constraints

  // bounds for decision variables
  double lb[] = {0.0, x1 * 0.95, u0}; /* lower bounds */
  double ub[] = {x1 * 0.9, x1, u1};   /* lower bounds */

  opt = nlopt_create(NLOPT_LN_COBYLA, n);
  // opt = nlopt_create(NLOPT_LD_LBFGS, n);
  // opt = nlopt_create(NLOPT_LD_SLSQP, n);

  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);

  nlopt_set_min_objective(opt, VdV, NULL);
  nlopt_set_xtol_rel(opt, 1e-6);
  // nlopt_add_inequality_constraint(opt, myconstraint, NULL, 1e-8);

  double minf;
  x[0] = ub[0];
  x[1] = lb[1];
  x[2] = u1 * 0.5;

  if (1)
  {
    nlopt_result status = nlopt_optimize(opt, x, &minf);
    if (status < 0)
    {
      printf("nlopt failed! Error code: %d\n", status);
      std::cout << "Optimization result: " << nlopt_result_to_string(status) << std::endl;
    }
    else
      printf("found minimum at %0.10g: x[3] = %.10g %.10g %.10g\n", minf, x[0], x[1], x[2]);
  }

  nlopt_destroy(opt);

  if (1) //create data file for LAMMPS
  {
    // interpolate
    int n = 500;
    size_t N = 10 * n;
    double result;

    // std::vector<double> lgrid, tgrid;
    double *lgrid = new double[N + 1];
    double *xgrid = new double[N + 1];

    double l1 = x[0], l2 = x[1], u2 = x[2];

    double a = 2 * (u2 - u0) / (l2 - l1) - (u1 - u2) / (x1 - l2);
    double b = ((u1 - u2) / (x1 - l2) - (u2 - u0) / (l2 - l1)) / (l2 - l1);

    for (int i = 0; i <= N; i++)
    {
      double xi = i * 1.0 * x1 / N;

      // tgrid.push_back(t); lgrid.push_back(result); c++;
      xgrid[i] = xi;
      double wi;
      if (xi < x[0])
      {
        result = xi;
        wi = u0;
      }
      else if (xi < x[1])
      {
        double t2 = 2 * b * (xi - l1);
        double L2 = ((a + t2) * sqrt(1 + (a + t2) * (a + t2)) + asinh(a + t2) - a * sqrt(1 + a * a) - asinh(a)) / (4 * b);

        result = l1 + L2;
        wi = u0 + a * (xi - l1) + b * (xi - l1) * (xi - l1);
      }
      else
      {
        double t2 = 2 * b * (l2 - l1);
        double L2 = ((a + t2) * sqrt(1 + (a + t2) * (a + t2)) + asinh(a + t2) - a * sqrt(1 + a * a) - asinh(a)) / (4 * b);

        wi = u2 + (u1 - u2) * (xi - l2) / (x1 - l2);
        double L3 = sqrt((wi - u2) * (wi - u2) + (xi - l2) * (xi - l2));
        result = l1 + L2 + L3;
      }

      lgrid[i] = result;

      printf("\n%lf %lf %lf", xgrid[i], lgrid[i], wi);
    }

    ///////////
    std::vector<Vector3> points;

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N + 1);
    gsl_spline_init(spline_steffen, lgrid, xgrid, N + 1);

    double dL = (lgrid[N] - lgrid[0]) / n;

    int Ni = n;
    // discrete.clear();
    for (int i = 0; i <= Ni; ++i)
    // for (float i = 0.5; i < Ni; ++i)
    {
      double li = (1 - i * 1.0 / Ni) * lgrid[0] + (i * 1.0 / Ni) * lgrid[N];
      double xi = gsl_spline_eval(spline_steffen, li, acc);

      double wi;
      if (xi < x[0])
        wi = u0;
      else if (xi < x[1])
        wi = u0 + a * (xi - l1) + b * (xi - l1) * (xi - l1);
      else
        wi = u2 + (u1 - u2) * (xi - l2) / (x1 - l2);

      printf("%g : %g %g\n", li, xi, wi);

      points.push_back(Vector3(xi, wi, 0.0));
    }

    for (int i = 0; i <= Ni; ++i)
    // for (float i = 0.5; i < Ni; ++i)
    {
      double li = (1 - i * 1.0 / Ni) * lgrid[0] + (i * 1.0 / Ni) * lgrid[N];
      double xi = gsl_spline_eval(spline_steffen, li, acc);

      double wi;
      if (xi < x[0])
        wi = u0;
      else if (xi < x[1])
        wi = u0 + a * (xi - l1) + b * (xi - l1) * (xi - l1);
      else
        wi = u2 + (u1 - u2) * (xi - l2) / (x1 - l2);

      //printf("%g : %g %g\n", li, xi, wi);

      points.push_back(Vector3(xi, -wi, 0.0));
    }


    gsl_spline_free(spline_steffen);
    gsl_interp_accel_free(acc);

    delete[] lgrid, xgrid;

    /////
    // Generate LAMMPS datafile

    FILE *fp1 = fopen("data.discrete", "w");
    int index;
    int nf = 2;
    double bounds = 1000;

    fprintf(fp1, "LAMMPS Description\n\n%ld atoms\n%d bonds\n%d angles\n0 dihedrals\n0 impropers\n\n", points.size(), 2*n, 2 * (n - 1));
    fprintf(fp1, "2 atom types\n1 bond types\n1 angle types\n\n");
    fprintf(fp1, "%lf %lf xlo xhi\n%lf %lf ylo yhi\n%lf %lf zlo zhi\n\n", -bounds, bounds, -bounds, bounds, -bounds, bounds);
    fprintf(fp1, "Masses\n\n1 1.0\n2 1.0\n\n");
    //fprintf(fp1, "Pair Coeffs\n\n1 %.2f %.2f %.2f\n\n", De * dL * dL, bm, re);

    fprintf(fp1, "Bond Coeffs\n\n1 %.4f %.6f\n\n", 0.5 * E * M_PI * R * R / dL, dL);
    fprintf(fp1, "Angle Coeffs\n\n1 %.4f 180.0\n\n", 0.25 * E * M_PI * R * R * R * R / dL);

    fprintf(fp1, "Atoms\n\n");
    index = 1;
    loop(i, points.size())
    {
      fprintf(fp1, "%d %d %d 0.00 %.6f %.6f %.6f\n", i + 1, i / (n + 1) + 1, i / (n + 1) + 1, points[i].comp[0], points[i].comp[1], points[i].comp[2]);
    }

    fprintf(fp1, "\nBonds\n\n");
    index = 1;
    loop(i, nf) loop2(j, 1, n + 1)
    {
      fprintf(fp1, "%d 1 %d %d\n", index++, i * (n + 1) + j, i * (n + 1) + j + 1);
    }

    fprintf(fp1, "\nAngles\n\n");
    index = 1;
    loop(i, nf) loop2(j, 1, n)
    {
      fprintf(fp1, "%d 1 %d %d %d\n", index++, i * (n + 1) + j, i * (n + 1) + j + 1, i * (n + 1) + j + 2);
    }

    fclose(fp1);
  }


  if (1) //create bezier fit
  {
    column_vector params(4);
    params = points1[NPOINTS-1].x, points1[NPOINTS-1].y, points1[0].x, points1[0].y;

    dlib::find_min_using_approximate_derivatives(
            dlib::cg_search_strategy(),
            dlib::objective_delta_stop_strategy(1).be_verbose(),
            rateCurve,
            params,
            -1);

    printf("p1x = %f;\n", points1[0].x);
    printf("p1y = %f;\n", points1[0].y);
    printf("c1x = %f;\n", params(0,0));
    printf("c1y = %f;\n", params(1,0));
    printf("c2x = %f;\n", params(2,0));
    printf("c2y = %f;\n", params(3,0));
    printf("p2x = %f;\n", points1[NPOINTS-1].x);
    printf("p2y = %f;\n", points1[NPOINTS-1].y);
  }

  delete[] x, df;
  return 0;
}

#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

//INCLUDE GSL HEADERS
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>

#include "cuba.h"
#include <nlopt.h>
#include "potential.h"

typedef std::stringstream sss;

#define ABS_ERR 0
#define REL_ERR 1e-4
#define REL_ERR2 1e-12
#define GRAD_TOL 1e-14
#define RTOL 1e-8

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

#define BTYPES 1
double EA[BTYPES] = {0}; //Axial elastic constants (EA)
double EI[BTYPES] = {0}; //Bending stiffness constants (EI)
double rho[BTYPES] = {0}; //mass density of bezier
double cfac = 0;

#include "util.h"
#include "cbezier.h"
#include "inter.h"
#include "fibnetwork.h"

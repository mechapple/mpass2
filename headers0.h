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

#include "cuba.h"
#include <nlopt.h>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

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

#define loop(x,n) for(int x = 0; x < n; ++x)
#define loop2(x,a,b) for(int x = a; x < b; ++x)

#define BTYPES 1
double EA[BTYPES] = {0}; //Axial elastic constants (EA)
double EI[BTYPES] = {0}; //Bending stiffness constants (EI)
double cfac = 50;

typedef std::stringstream sss;

struct Vector3 {
  double comp[3];
  
  Vector3(double x, double y, double z)
  {
    comp[0] = x;  comp[1] = y;  comp[2] = z;
  }
  
  Vector3()
  {
    comp[0] = 0;  comp[1] = 0;  comp[2] = 0;
  }
};

struct iVector3 {
  int comp[3];
  
  iVector3(int x, int y, int z)
  {
    comp[0] = x;  comp[1] = y;  comp[2] = z;
  }
  
  iVector3()
  {
    comp[0] = 0;  comp[1] = 0;  comp[2] = 0;
  }
};

int gsl_matrix_print(const gsl_matrix *m, const char *title)
{
  int status, n = 0;
  printf("%s\n",title);
  
  for (size_t i = 0; i < m->size1; i++) {
    for (size_t j = 0; j < m->size2; j++) {
      if ((status = printf("%g ", gsl_matrix_get(m, i, j))) < 0)
        return -1;
      n += status;
    }
  
    if ((status = printf("\n")) < 0)
      return -1;
    n += status;
  }
  
  return n;
}

class cBezier {
  
  public:
    double x[12],v[12],f[12];
    int CP[4], type, id;
    double energy, length, length0, axialE, bendE, cohE;
    int leftB, rightB;
    
    gsl_matrix * P = gsl_matrix_alloc (4, 3);
    gsl_matrix * PB = gsl_matrix_alloc (3, 4);
    gsl_matrix * P_B = gsl_matrix_alloc (4, 4);
		gsl_matrix * PP = gsl_matrix_alloc (4, 4);
  
    void interpolate(std::vector<Vector3> &, int);
    void initialize_bezier();
    void update_bezier();
    void axial_ef();
    void bending_ef();
    
  private:
    void length_bezier();
    void set_length0();
    
};

void cBezier::interpolate(std::vector<Vector3> &points, int n) {
  
  for(double t=1.0/n; t<(1.0+1.0/n); t+=1.0/n)
  {
    double tc = 1.0-t;
    double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
    
    double rx = X[0]*x[0] + X[1]*x[3] + X[2]*x[6] + X[3]*x[9];
    double ry = X[0]*x[1] + X[1]*x[4] + X[2]*x[7] + X[3]*x[10];
    double rz = X[0]*x[2] + X[1]*x[5] + X[2]*x[8] + X[3]*x[11];
    
    points.push_back(Vector3(rx,ry,rz));
  }
}

void cBezier::initialize_bezier() 
{
  double Bvec[] = {-1, 3, -3, 1, 3, -6, 3, 0, -3, 3, 0, 0, 1, 0, 0, 0};
  
  const gsl_matrix_const_view P1 = gsl_matrix_const_view_array( x, 4, 3 );
	const gsl_matrix_const_view B1 = gsl_matrix_const_view_array( Bvec, 4, 4 );
  
  gsl_matrix_memcpy(P, &P1.matrix);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, &P1.matrix, &B1.matrix, 0.0, PB);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, PB, PB, 0.0, P_B);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, &P1.matrix, &P1.matrix, 0.0, PP);
  
  gsl_matrix_print(PP,"PP:");
  length_bezier(); set_length0();
  printf("Length = %lf\n",length);
  
  loop(i,12) {v[i] = 0.0, f[i] = 0.0;}
}

void cBezier::update_bezier() 
{
  double Bvec[] = {-1, 3, -3, 1, 3, -6, 3, 0, -3, 3, 0, 0, 1, 0, 0, 0};
  
  const gsl_matrix_const_view P1 = gsl_matrix_const_view_array( x, 4, 3 );
	const gsl_matrix_const_view B1 = gsl_matrix_const_view_array( Bvec, 4, 4 );
  
  gsl_matrix_memcpy(P, &P1.matrix);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, &P1.matrix, &B1.matrix, 0.0, PB);
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, PB, PB, 0.0, P_B);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, &P1.matrix, &P1.matrix, 0.0, PP);
  
  //gsl_matrix_print(PP,"PP:");
  //length_bezier(); set_length0();
  //printf("Length = %lf\n",length);
  
  loop(i,12) {v[i] = 0.0, f[i] = 0.0;}
}

double ds(double t, void *p)
{
	cBezier params = *(cBezier *) p;
	
  double tc = 1.0-t;
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
	double J = 0;
	
  loop(i,4)
		loop(j,4)
			J += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
			
	return sqrt(J);	
}

void cBezier::length_bezier()
{
  double result, error; size_t nev;
  gsl_function B1;
  B1.function = &ds;
  B1.params = this;
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(100);
  
  gsl_integration_cquad (&B1, 0.0, 1.0, ABS_ERR, REL_ERR2,w1, &result, &error, &nev);
  length = result;
}

void cBezier::set_length0() 
{
  length0 = length;
}

static int dJ(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0;
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  
  loop(i,4) loop(j,3) {
    double sum = 0.0;
    loop(k,4) sum += DX[i]*DX[k]*gsl_matrix_get(params.P, k, j);
    ff[i*3+j] = sum/sqrt(Jsq);
  }
  
  return 0;
}

void cBezier::axial_ef()
{
  length_bezier();
  axialE = 0.5*EA[type-1]*(length-length0)*(length-length0);
  
  int comp, nregions, neval, fail;
  double integral[12], error[12], prob[12];
  
  Cuhre(2, 12, dJ, this, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
    
  loop(i,12) {
    f[i] += -0.5*EA[type-1]*(length-length0)*integral[i];
    //f[i] += integral[i];
  }
}

double kappa_squared_ds(double t, void *p)
{
	cBezier params = *(cBezier *) p;
	
  double tc = 1.0-t;
  
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double Y[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  double Z[4] = {6*tc, 3*(2-6*tc), 3*(2-6*t), 6*t};
    
  double dr2=0, ddr2=0, drddr=0;
  loop(i,4) loop(j,4) {
    dr2 += Y[i]*gsl_matrix_get(params.PP, i, j)*Y[j];
    ddr2 += Z[i]*gsl_matrix_get(params.PP, i, j)*Z[j];
    drddr += Y[i]*gsl_matrix_get(params.PP, i, j)*Z[j];
  }
  
  double drxddr2 = dr2*ddr2-drddr*drddr;
  double kappa2 = drxddr2/(dr2*dr2*dr2);
  
	return kappa2*sqrt(dr2);	
}

static int dB(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0;
  
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double Y[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  double Z[4] = {6*tc, 3*(2-6*tc), 3*(2-6*t), 6*t};
    
  double dr2=0, ddr2=0, drddr=0;
  loop(i,4) loop(j,4) {
    dr2 += Y[i]*gsl_matrix_get(params.PP, i, j)*Y[j];
    ddr2 += Z[i]*gsl_matrix_get(params.PP, i, j)*Z[j];
    drddr += Y[i]*gsl_matrix_get(params.PP, i, j)*Z[j];
  }
  
  double drxddr2 = dr2*ddr2-drddr*drddr;
  double kappa2 = drxddr2/(dr2*dr2*dr2);
  
  loop(i,4) loop(j,3) {
    double sum1 = 0.0, sum2 = 0.0;
    loop(k,4) { 
      sum1 += 2*dr2*dr2*dr2*(  dr2*Z[i]*Z[k] + ddr2*Y[i]*Y[k] - drddr*(Y[i]*Z[k]+Z[i]*Y[k]) )*gsl_matrix_get(params.P, k, j);
      sum1 -= 6*drxddr2*dr2*dr2*Y[i]*Y[k]*gsl_matrix_get(params.P, k, j);
      sum2 += Y[i]*Y[k]*gsl_matrix_get(params.P, k, j);
    }
    
    ff[i*3+j] = (sum1*dr2/pow(dr2,6) + sum2*kappa2)/sqrt(dr2);
  }
  
  return 0;
}

void cBezier::bending_ef()
{
  double resultc, errorc; size_t nev;
  gsl_function B1;
  B1.function = &kappa_squared_ds;
  B1.params = this;
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(100);
  
  gsl_integration_cquad (&B1, 0.0, 1.0, ABS_ERR, REL_ERR2,w1, &resultc, &errorc, &nev);
  bendE = 0.5*EI[type-1]*resultc;

  int comp, nregions, neval, fail;
  double integral[12], error[12], prob[12];
  
  Cuhre(2, 12, dB, this, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
    
  loop(i,12) {
    f[i] += -0.5*EI[type-1]*integral[i];
  }
}

class fibnetwork {
  
  public:
    int np,nb,nf,nt,ncc;
    double bounds[3][2];
    double axialE, bendE, cbendE;
    
    std::map<int,double> fixed_dof;
    std::vector<iVector3> C1_cps;
    std::vector<Vector3> CPx,CPv,CPf;
    std::vector<cBezier> Bz_list;
    std::vector<std::vector<int>> fibers; 

    void readfile(char *);
    void printdata();
    void printlammps(char *,char *, int);
    void printlammps_cps(char *,char *);
    void compute_ef();
    void constraint_energy();
    void minimize();
};

void fibnetwork::readfile(char *filename)
{
  std::ifstream infile(filename);
  
  std::string word,line;
  loop(i,2) std::getline(infile, line);
  
  sss ss0;
  std::getline(infile, line); ss0 << line; ss0 >> np; ss0.str("");
  std::getline(infile, line); ss0 << line; ss0 >> nb; ss0.str("");
  std::getline(infile, line); ss0 << line; ss0 >> nf; ss0.str("");
  std::getline(infile, line); ss0 << line; ss0 >> ncc; ss0.str("");
  std::getline(infile, line);
  std::getline(infile, line); ss0 << line; ss0 >> nt; ss0.str("");
  std::getline(infile, line);
  
  loop(i,3) {
    sss ss;
    std::getline(infile, line); ss << line; 
    ss >> bounds[i][0] >> bounds[i][1];    
  }
  loop(i,7) std::getline(infile, line);
  
  loop(i,BTYPES) {
    sss ss;
    std::getline(infile, line); ss << line;
    int ctmp;  ss >> ctmp >> EA[ctmp-1] >> EI[ctmp-1]; ss.str("");
    //printf("%s %d %lf %lf\n",line.c_str(),ctmp,EA[ctmp-1],EI[ctmp-1]);
  }
  
  loop(i,3) std::getline(infile, line);
  
  loop(i,np) {
    sss ss;
    double x[3]; int tmp,fdof[3];
    std::getline(infile, line); ss << line; 
    ss >> tmp >> x[0] >> x[1] >> x[2] >> fdof[0] >> fdof[1] >> fdof[2];
    
    loop(j,3) if(fdof[j]==0) fixed_dof.insert(std::make_pair((tmp-1)*3+j+1, x[j] ));
    
    double rn = RTOL*((double) std::rand() / (RAND_MAX));
    
    CPx.push_back(Vector3(x[0]+rn,x[1]+rn,x[2]+rn));
    CPv.push_back(Vector3(0,0,0));
    CPf.push_back(Vector3(0,0,0));
  }
  loop(i,3) std::getline(infile, line);
  
  loop(i,nb) {
    sss ss;
    cBezier bzi; int tmp[6]; 
    std::getline(infile, line); ss << line; 
    ss >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3] >> tmp[4] >> tmp[5];
    
    bzi.type = tmp[1]; bzi.id = tmp[0];
    loop(j,4) {
      bzi.CP[j] = tmp[j+2]; 
      int cpid = bzi.CP[j]-1;
      loop(k,3) bzi.x[k+j*3] = CPx[cpid].comp[k];
    }
    bzi.initialize_bezier();
    
    Bz_list.push_back(bzi);
  }
  loop(i,3) std::getline(infile, line);
  
  loop(i,nf) {
    sss ss; int tmp;
    std::vector<int> fbz;
    std::getline(infile, line); ss << line; ss >> tmp;
    while (ss >> tmp) fbz.push_back(tmp);
    
    loop(j,fbz.size()) {
      if(j==0) {Bz_list[fbz[j]-1].leftB = -1; Bz_list[fbz[j]-1].rightB = fbz[j+1];}
      else if(j== fbz.size()-1) {Bz_list[fbz[j]-1].leftB = fbz[j-1]; Bz_list[fbz[j]-1].rightB = -1;}
      else {Bz_list[fbz[j]-1].leftB = fbz[j-1]; Bz_list[fbz[j]-1].rightB = fbz[j+1];}
    } 
    
    fibers.push_back(fbz);
  }
  loop(i,3) std::getline(infile, line);
  
  loop(i,ncc) {
    sss ss;
    cBezier bzi; int tmp[4]; 
    std::getline(infile, line); ss << line; 
    ss >> tmp[0] >> tmp[1] >> tmp[2] >> tmp[3];
    
    C1_cps.push_back(iVector3(tmp[1],tmp[2],tmp[3]));
  }
  
  CPv.resize(CPx.size()); CPf.resize(CPx.size());
  
  infile.close();
}

void fibnetwork::printdata()
{
  printf("Control points \n");
  loop(i,np) {
    printf("%d %lf %lf %lf \n",i+1, CPx[i].comp[0], CPx[i].comp[1], CPx[i].comp[2]);
  }
  
  printf("\nCubic Beziers \n");
  loop(i,nb) {
    printf("%d %d %d %d %d %d (%d,%d)\n",Bz_list[i].id, Bz_list[i].type, Bz_list[i].CP[0], Bz_list[i].CP[1], Bz_list[i].CP[2], Bz_list[i].CP[3], Bz_list[i].leftB, Bz_list[i].rightB);
  }
  
  printf("\nFibers \n");
  loop(i,nf) {
    std::cout << i+1 << ' ';
    for (auto k: fibers[i]) std::cout << k << ' ';
    std::cout << std::endl ;
  }
}

void fibnetwork::printlammps_cps(char *filename, char *mode)
{
  FILE *fp = fopen(filename, mode);
  
  fprintf(fp, "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n%d\n",np);
  fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n%2.6f %2.6f \n%2.6f %2.6f \n%2.6f %2.6f \n", 
      bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1], bounds[2][0], bounds[2][1]);
  fprintf(fp, "ITEM: ATOMS id mol type q xu yu zu\n");
  
  loop(i,np) {
    fprintf(fp, "%d 1 1 0.00 %.6f %.6f %.6f\n", i+1, CPx[i].comp[0], CPx[i].comp[1], CPx[i].comp[2]);
  }
  
  fclose(fp);
}

void fibnetwork::printlammps(char *filename, char *mode, int n)
{
  //int num_coord = 3*(n*nb+1);
  
  std::vector<Vector3> points;
  points.push_back(Vector3(Bz_list[fibers[0][0]-1].x[0],Bz_list[fibers[0][0]-1].x[1],Bz_list[fibers[0][0]-1].x[2]));
  
  loop(i,nf) {
    for (auto k: fibers[i])
    {
      Bz_list[k-1].interpolate(points,n);
    }
  }  
  
  FILE *fp = fopen(filename, mode);
  
  fprintf(fp, "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n%ld\n",points.size());
  fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n%2.6f %2.6f \n%2.6f %2.6f \n%2.6f %2.6f \n", 
      bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1], bounds[2][0], bounds[2][1]);
  fprintf(fp, "ITEM: ATOMS id mol type q xu yu zu\n");
  
  loop(i,points.size()) {
    fprintf(fp, "%d 1 1 0.00 %.6f %.6f %.6f\n", i+1, points[i].comp[0], points[i].comp[1], points[i].comp[2]);
  }
  
  fclose(fp);
}

void fibnetwork::constraint_energy()
{
  cbendE = 0.0;
  
  loop(i,C1_cps.size()) {
    int a = C1_cps[i].comp[0];
    int b = C1_cps[i].comp[1];
    int c = C1_cps[i].comp[2];
    
    Vector3 p,q;
    
    loop(j,3) {
      p.comp[j] = CPx[a-1].comp[j] - CPx[b-1].comp[j];
      q.comp[j] = CPx[c-1].comp[j] - CPx[b-1].comp[j];
    }

    double pp=0,qq=0,pq=0;
    loop(j,3) {
      pp += p.comp[j]*p.comp[j];
      qq += q.comp[j]*q.comp[j];
      pq += p.comp[j]*q.comp[j];
    }
    
    cbendE += 0.5*cfac*EI[0]*(pp*qq - pq*pq); 
    
    if(1) {
      loop(j,3) {
        CPf[a-1].comp[j] -= EI[0]*cfac*(2*qq*p.comp[j] - 2*pq*q.comp[j]);
        CPf[c-1].comp[j] -= EI[0]*cfac*(2*pp*q.comp[j] - 2*pq*p.comp[j]);
        CPf[b-1].comp[j] -= EI[0]*cfac*(2*pq*(p.comp[j]+q.comp[j]) - 2*pp*q.comp[j] - 2*qq*p.comp[j]);   
      }
    }
  }
}

void fibnetwork::compute_ef()
{
  axialE = 0.0; bendE = 0.0;
  loop(i,np) loop(j,3) CPf[i].comp[j] = 0.0;
  
  loop(i,nf) {
    for (auto k: fibers[i])
    {
      loop(b,4) loop(j,3) Bz_list[k-1].f[b*3+j] = 0.0;
      
      Bz_list[k-1].axial_ef();      
      //printf("Forces: \n"); loop(i,4) printf("%.4e %.4e %.4e \n", Bz_list[k-1].f[i*3],Bz_list[k-1].f[i*3+1],Bz_list[k-1].f[i*3+2]);
      Bz_list[k-1].bending_ef();
      //printf("Forces: \n"); loop(i,4) printf("%.4e %.4e %.4e \n", Bz_list[k-1].f[i*3],Bz_list[k-1].f[i*3+1],Bz_list[k-1].f[i*3+2]);
      //printf("%d %d %lf %lf\n",i+1,k,Bz_list[k-1].axialE, Bz_list[k-1].bendE);
      
      loop(b,4) {
        int id = Bz_list[k-1].CP[b];
        loop(j,3) CPf[id-1].comp[j] += Bz_list[k-1].f[b*3+j];
      }
      
      axialE += Bz_list[k-1].axialE;
      bendE += Bz_list[k-1].bendE;
    }
  }
  constraint_energy();  
  printf("Energies: %lf %lf %lf\n", axialE, bendE, cbendE);
  
}

double VdV(unsigned n, const double *input, double *grad, void *fn_data)
{
  fibnetwork *FN = (fibnetwork *) fn_data;
    
  loop(i,FN->np) loop(j,3) FN->CPx[i].comp[j] = input[i*3+j];
    
  loop(i,FN->nb) {
    
    loop(j,4) {
      
      int cpid = FN->Bz_list[i].CP[j]-1;
      loop(k,3) 
        FN->Bz_list[i].x[j*3+k] = FN->CPx[cpid].comp[k];
    }
    FN->Bz_list[i].update_bezier();
  }
  
  FN->compute_ef();
      
  if (grad) {
    loop(i,FN->np) loop(j,3) grad[i*3+j] = -FN->CPf[i].comp[j] ;
  }
  
  FN->printlammps_cps((char*) "Output_cps.lammpstrj",(char*) "a");
  FN->printlammps((char*) "Output.lammpstrj",(char*) "a",20);
  
  return (FN->axialE + FN->bendE);
}

void fixedconstraints(unsigned m, double *result, unsigned n, const double* input, double* grad, void* fn_data)
{
  fibnetwork *FN = (fibnetwork *) fn_data;
  
  if (grad) loop(i,m*n) grad[i] = 0.0;
  
  int ci = 0;
  for(std::map<int,double>::iterator it = FN->fixed_dof.begin(); it != FN->fixed_dof.end(); it++) {
    result[ci] = input[it->first-1] - it->second;
    if (grad) grad[ci*n+(it->first-1)] = 1.0;
    ci++;
  }
}

void C1_continuity(unsigned m, double *result, unsigned n, const double* input, double* grad, void* fn_data)
{
  fibnetwork *FN = (fibnetwork *) fn_data;
  
  if (grad) loop(i,m*n) grad[i] = 0.0;
  
  loop(i,FN->C1_cps.size()) {
    int a = FN->C1_cps[i].comp[0];
    int b = FN->C1_cps[i].comp[1];
    int c = FN->C1_cps[i].comp[2];
    
    Vector3 p,q;
    
    loop(j,3) {
      p.comp[j] = input[(a-1)*3+j] - input[(b-1)*3+j];
      q.comp[j] = input[(c-1)*3+j] - input[(b-1)*3+j];
    }

    double pp=0,qq=0,pq=0;
    loop(j,3) {
      pp += p.comp[j]*p.comp[j];
      qq += q.comp[j]*q.comp[j];
      pq += p.comp[j]*q.comp[j];
    }
    
    result[i] = pp*qq - pq*pq; 
    
    if(grad) {
      loop(j,3) {
        grad[i*n + (a-1)*3+j] = 2*qq*p.comp[j] - 2*pq*q.comp[j];
        grad[i*n + (c-1)*3+j] = 2*pp*q.comp[j] - 2*pq*p.comp[j];
        grad[i*n + (b-1)*3+j] = 2*pq*(p.comp[j]+q.comp[j]) - 2*pp*q.comp[j] - 2*qq*p.comp[j];    
        
      }
    }
  }
}

void fibnetwork::minimize()
{
  double *x = new double[np*3];
  double *df = new double[np*3];
  loop(i,np) loop(j,3) x[i*3+j] = CPx[i].comp[j];
  
  nlopt_opt opt;

  //establish sizes
  unsigned n = np*3; //number of decision variables
  unsigned m_in = 0; //number of inequality constraints

  //bounds for decision variables
  //double lb[] = { 0.3, -HUGE_VAL, -HUGE_VAL,  -HUGE_VAL, -HUGE_VAL }; /* lower bounds */
  //double ub[] = { HUGE_VAL, HUGE_VAL, 5, HUGE_VAL, HUGE_VAL }; /* lower bounds */
  
  double *lb = new double[np*3];
  double *ub = new double[np*3];
  
  loop(i,np) loop(j,3) {
    lb[i*3+j] = bounds[j][0];
    ub[i*3+j] = bounds[j][1];
  }
  
  opt = nlopt_create(NLOPT_LN_COBYLA, n);
  //opt = nlopt_create(NLOPT_LD_LBFGS, n);
  //opt = nlopt_create(NLOPT_LD_SLSQP, n);
  
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);

  nlopt_set_min_objective(opt, VdV, this);
  
  //double tol_eq[]={1e-8,1e-8};
  unsigned m_eq = fixed_dof.size(); //number of equality constraints
  double *fctol = new double[m_eq];  loop(i,m_eq) fctol[i] = 1e-8;
  nlopt_add_equality_mconstraint(opt, m_eq, fixedconstraints, this, fctol);

  nlopt_set_xtol_rel(opt, 1e-4);
  
  double minf;
  
  if(1) {
    nlopt_result status = nlopt_optimize(opt, x, &minf);
    if (status < 0) {
      printf("nlopt failed! Error code: %d\n",status);
      std::cout << "Optimization result: " << nlopt_result_to_string(status) << std::endl;
    }
    else
      printf("found minimum at %0.10g\n", minf);  
  }
  
  unsigned m_eq2 = C1_cps.size(); //number of equality constraints
  double *cctol = new double[m_eq2];  loop(i,m_eq2) cctol[i] = 1e-8;
  nlopt_add_equality_mconstraint(opt, m_eq2, C1_continuity, this, cctol);
  
  cfac = 0.0;
  
  if(1) {
    nlopt_result status = nlopt_optimize(opt, x, &minf);
    if (status < 0) {
      printf("nlopt failed! Error code: %d\n",status);
      std::cout << "Optimization result: " << nlopt_result_to_string(status) << std::endl;
    }
    else
      printf("found minimum at %0.10g\n", minf);  
  }
  
  nlopt_destroy(opt);
  //double fn_energy = VdV(np, x, df, this);
  //printf("FN Energy: %lf\nFN Gradients: \n",fn_energy);
  //loop(i,np) printf("%.4e %.4e %.4e\n",df[i*3],df[i*3+1],df[i*3+2]);
  
  //printf("Fixed dofs: \n");
  //for(std::map<int,double>::iterator it = fixed_dof.begin(); it != fixed_dof.end(); it++)
    //printf("%d %lf\n",it->first,it->second);
    
  //printf("Collinear cps: \n");
  //loop(i,C1_cps.size()) {
    //int a = C1_cps[i].comp[0];
    //int b = C1_cps[i].comp[1];
    //int c = C1_cps[i].comp[2];
    //printf("%d %d %d \n",a,b,c);
  //}
  
  delete [] lb,ub,x,df,fctol;
}

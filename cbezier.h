class cBezier {
  
  public:
    double x[12],v[12],f[12],a[12];
    int CP[4], type, id, fiber;
    double energy, length, length0, axialE, bendE, cohE, kinE;
    int leftB, rightB;
    
    gsl_matrix * P = gsl_matrix_alloc (4, 3);
    gsl_matrix * PB = gsl_matrix_alloc (3, 4);
    gsl_matrix * P_B = gsl_matrix_alloc (4, 4);
		gsl_matrix * PP = gsl_matrix_alloc (4, 4);
    
    gsl_matrix * Q = gsl_matrix_alloc (4, 3);
    gsl_matrix * PQ_symm = gsl_matrix_alloc (4, 4);
    gsl_matrix * QQ = gsl_matrix_alloc (4, 4);
    gsl_matrix * M = gsl_matrix_alloc (4, 4);
    gsl_matrix * Mdot = gsl_matrix_alloc (4, 4);
    gsl_matrix * N = gsl_matrix_alloc (4, 4);
  
    void interpolateT(std::vector<Vector3> &, int);
    void interpolateL(std::vector<Vector3> &, int);
    void initialize_bezier();
    void update_bezier();
    void axial_ef();
    void bending_ef();
    void kin_energy();
    
  private:
    void length_bezier();
    void set_length0();
    
};

void cBezier::interpolateT(std::vector<Vector3> &points, int n) {
  
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

void cBezier::interpolateL(std::vector<Vector3> &points, int n) {
  
  double result, error; size_t nev;
  
  gsl_function B1;  B1.function = &ds;  B1.params = this;
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(100);

  //std::map<double, double> TLmap;
  size_t N = 10*n+1;
  double *lgrid = new double[N];
  double *tgrid = new double[N];
  
  int c=0;
  
  for(double t=0.0; t<(1.0+0.1/n); t+=0.1/n)
  {
    gsl_integration_cquad (&B1, 0.0, t, ABS_ERR, REL_ERR2,w1, &result, &error, &nev);
    //TLmap.insert(std::make_pair(t, result));
    
    tgrid[c] = t; lgrid[c] = result;  c++;
  }
  
  gsl_integration_cquad_workspace_free (w1);
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N);
  gsl_spline_init(spline_steffen, lgrid, tgrid, N);

  int Ni = n;
  for (int i = 1; i <= Ni; ++i)
  {
    double li = (1 - i*1.0 / Ni) * lgrid[0] + (i*1.0 / Ni) * lgrid[N-1];
    double t = gsl_spline_eval(spline_steffen, li, acc);

    //printf("%g : %g\n", li, t);
    
    double tc = 1.0-t;
    double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
    
    double rx = X[0]*x[0] + X[1]*x[3] + X[2]*x[6] + X[3]*x[9];
    double ry = X[0]*x[1] + X[1]*x[4] + X[2]*x[7] + X[3]*x[10];
    double rz = X[0]*x[2] + X[1]*x[5] + X[2]*x[8] + X[3]*x[11];
    
    points.push_back(Vector3(rx,ry,rz));
  }

  gsl_spline_free(spline_steffen);
  gsl_interp_accel_free(acc);

  delete lgrid,tgrid;
  
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
  
  //gsl_matrix_print(PP,"PP:");
  length_bezier(); set_length0();
  //printf("Length = %lf\n",length);
  
  loop(i,12) {v[i] = 0.0, f[i] = 0.0;}
  
  update_bezier();
}

static int XXJ(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0;
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  double J = sqrt(Jsq);
  
  loop(i,4) loop(j,4) {
    ff[i*4+j] = X[i]*X[j]*J;
  }
  
  return 0;
}

static int XXf(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0, f = 0;
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  double J = sqrt(Jsq);
  
  loop(i,4) loop(j,4) f += DX[i]*gsl_matrix_get(params.PQ_symm, i, j)*DX[j];
  f = f/2;
  
  loop(i,4) loop(j,4) {
    ff[i*4+j] = X[i]*X[j]*f/J;
  }
  
  return 0;
}

static int YYf(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *p) {
  
  cBezier params = *(cBezier *) p;

  double t = xx[0], tc = 1.0 - xx[0], Jsq = 0, f = 0;
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  loop(i,4) loop(j,4) Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
  double J = sqrt(Jsq);
  
  loop(i,4) loop(j,4) f += X[i]*gsl_matrix_get(params.QQ, i, j)*X[j];
  f = f/2;
  
  loop(i,4) loop(j,4) {
    ff[i*4+j] = DX[i]*DX[j]*f/J;
  }
  
  return 0;
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
  
  //update M and Mdot
  int comp, nregions, neval, fail;
  double integral[16], error[16], prob[16];
  
  Cuhre(2, 16, XXJ, this, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
  
  loop(i,4) loop(j,4) gsl_matrix_set(M,i,j,rho[type-1]*(length0/length)*integral[i*4+j]);
  
  const gsl_matrix_const_view Q1 = gsl_matrix_const_view_array( v, 4, 3 );
  gsl_matrix_memcpy(Q, &Q1.matrix);
  
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Q, P, 0.0, PQ_symm);
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, P, Q, 1.0, PQ_symm);
    
  Cuhre(2, 16, XXf, this, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
    
  loop(i,4) loop(j,4) gsl_matrix_set(Mdot,i,j,rho[type-1]*(length0/length)*integral[i*4+j]);  
  
  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Q, Q, 0.0, QQ);
  
  Cuhre(2, 16, YYf, this, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
    
  loop(i,4) loop(j,4) gsl_matrix_set(N,i,j,rho[type-1]*(length0/length)*integral[i*4+j]);  
  
  //gsl_matrix_print(PP,"PP:");
  //length_bezier(); set_length0();
  //printf("Length = %lf\n",length);
  
  loop(i,12) {v[i] = 0.0, f[i] = 0.0;}
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

double vel_squared_ds(double t, void *p)
{
	cBezier params = *(cBezier *) p;
	
  double tc = 1.0-t, Jsq = 0;
    
  double X[4] = {tc*tc*tc, 3*tc*tc*t, 3*tc*t*t, t*t*t};
  double DX[4] = {-3*tc*tc, 3*tc*(tc-2*t), -3*t*(t-2*tc), 3*t*t};
  
  double vv=0;
  loop(i,4) loop(j,4) {
    Jsq += DX[i]*gsl_matrix_get(params.PP, i, j)*DX[j];
    vv += X[i]*gsl_matrix_get(params.QQ, i, j)*X[j];
  }
  
	return vv*sqrt(Jsq);	
}

void cBezier::kin_energy()
{
  double resultc, errorc; size_t nev;
  
  gsl_function B1;
  B1.function = &vel_squared_ds;
  B1.params = this;
  
  gsl_integration_cquad_workspace * w1 = gsl_integration_cquad_workspace_alloc(100);
  gsl_integration_cquad (&B1, 0.0, 1.0, ABS_ERR, REL_ERR2,w1, &resultc, &errorc, &nev);
  
  kinE = 0.5*rho[type-1]*(length0/length)*resultc;
}

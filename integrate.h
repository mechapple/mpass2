int func(double t, const double W[], double f[], void *fn_data)
{
  (void)(t); /* avoid unused parameter warning */
  fibnetwork *FN = (fibnetwork *) fn_data;

  int np = FN->np;	
  int nb = FN->nb;
  int nf = FN->nf;
  
  loop(i,np) loop(j,3) FN->CPx[i].comp[j] = W[i*3+j];
  loop(i,np) loop(j,3) FN->CPv[i].comp[j] = W[i*3+j+np*3];
  
  //apply fixed constraints
  for(std::map<int,double>::iterator it = FN->fixed_dof.begin(); it != FN->fixed_dof.end(); it++) {
    int index = it->first-1;
    int i0 = index/3;
    int j0 = index%3;
    FN->CPx[i0].comp[j0] = it->second;
  }
  
  if(1)
    loop(i,FN->C2_cps.size()) {
      int a = FN->C2_cps[i].comp[0];
      int b = FN->C2_cps[i].comp[1];
      int c = FN->C2_cps[i].comp[2];
      int d = FN->C2_cps[i].comp[3];
      int e = FN->C2_cps[i].comp[4];
      
      Vector3 db(FN->CPx[d-1].comp[0]-FN->CPx[b-1].comp[0], FN->CPx[d-1].comp[1]-FN->CPx[b-1].comp[1] , FN->CPx[d-1].comp[2]-FN->CPx[b-1].comp[2] );
      Vector3 eb(FN->CPx[e-1].comp[0]-FN->CPx[b-1].comp[0], FN->CPx[e-1].comp[1]-FN->CPx[b-1].comp[1] , FN->CPx[e-1].comp[2]-FN->CPx[b-1].comp[2] );
      Vector3 ad(FN->CPx[a-1].comp[0]-FN->CPx[d-1].comp[0], FN->CPx[a-1].comp[1]-FN->CPx[d-1].comp[1] , FN->CPx[a-1].comp[2]-FN->CPx[d-1].comp[2] );
      
      double m = sqrt(cross_mag(db,eb)/cross_mag(db,ad));
      double x = 1.0/(1.0+m);
      
      //EDIT to avoid singularities
      if(x<0.05) x = 0.05; if(x>0.95) x = 0.95;
      
      //x = 0.5;
      
      //printf("%d %d %d %d %d %d %lf %lf\n",i,a,b,c,d,e,m,x);
      loop(j,3) FN->CPx[c-1].comp[j] = x*FN->CPx[d-1].comp[j] + (1.0-x)*FN->CPx[b-1].comp[j];
    }
      
  loop(i,nb) {
    loop(j,4) {
      int cpid = FN->Bz_list[i].CP[j]-1;
      loop(k,3) {
        FN->Bz_list[i].x[j*3+k] = FN->CPx[cpid].comp[k];
        FN->Bz_list[i].v[j*3+k] = FN->CPv[cpid].comp[k];
	  }
    }
    FN->Bz_list[i].update_bezier();
  }
  
  loop(i,np) loop(j,3) {
    gsl_matrix_set(FN->Psys, i,j, FN->CPx[i].comp[j]);
    gsl_matrix_set(FN->Qsys, i,j, FN->CPv[i].comp[j]);
  }
  
  gsl_matrix *K1 = gsl_matrix_alloc (2*np, 3);
  
  FN->computeK(K1);

  loop(i,2*np) loop(j,3) f[i*3+j] = gsl_matrix_get(K1, i,j);
  
  gsl_matrix_free(K1);
  
  return GSL_SUCCESS;
}

void integrate(fibnetwork &fn, int nsteps, double tstep)
{
	size_t np = fn.np;
	
	fn.Psys = gsl_matrix_alloc (np, 3);
	fn.Qsys = gsl_matrix_alloc (np, 3);
	fn.Rsys = gsl_matrix_alloc (np, 3);
	fn.Fsys = gsl_matrix_alloc (np, 3);
	
	fn.GM = gsl_matrix_alloc (np, np);
	fn.GMdot = gsl_matrix_alloc (np, np);
	fn.GN = gsl_matrix_alloc (np, np);
  
	double *W = new double[2*np*3];
	loop(i,np) loop(j,3) W[i*3+j] = fn.CPx[i].comp[j];
	loop(i,np) loop(j,3) W[i*3+j+np*3] = fn.CPv[i].comp[j];
	
	gsl_odeiv2_system sys = {func, NULL, 2*np*3, &fn};

	gsl_odeiv2_driver * d = 
	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
	
	double t = 0.0, t1 = 100.0;
	
	for (int s = 1; s <= nsteps; s++)
    {
		double ti = tstep * s;
		
		int status = gsl_odeiv2_driver_apply (d, &t, ti, W);
		
		if (status != GSL_SUCCESS)
        {
			printf ("error, return value=%d\n", status);
			break;
        }
        
        //loop(i,np) loop(j,3) std::cout << " " << fn.CPx[i].comp[j];
        //std::cout << s << std::endl;
        
      fn.printlammps_cps((char*) "Output_cps.lammpstrj",(char*) "a");
      fn.printlammps((char*) "Output.lammpstrj",(char*) "a",20);
		
      printf("%d %lf %lf %lf %lf %lf\n", s, fn.axialE, fn.bendE, fn.cohE, fn.kinE, fn.totE);
	}
	
	//loop(i,2*np*3) W[i] += 0.001;
	
	loop(i,np) loop(j,3) fn.CPx[i].comp[j] = W[i*3+j];
	loop(i,np) loop(j,3) fn.CPv[i].comp[j] = W[i*3+j+np*3];
	
}

#define LJ_EPS 0.50
#define LJ_SIG 0.008

double V(double r) {
  double fsr = fabs(LJ_SIG/r);
	double rn = pow(fsr,12);
	return 4*LJ_EPS*(rn-fsr);
	//return 1.0;	
}

//f(r) = V'(r)/r
double f(double r) {
	
  double fsr = fabs(LJ_SIG/r);
	double rn = pow(fsr,12);
	
  return 4*LJ_EPS*(-12*rn+fsr)/(r*r);	
	//return 1.0;	
}

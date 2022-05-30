#define LJ_EPS 5.0
#define LJ_SIG 0.008

double V(double r) {
	double rn = pow(LJ_SIG/r,6);
	return 4*LJ_EPS*(rn*rn-rn);
	//return 1.0;	
}

//f(r) = V'(r)/r
double f(double r) {
	double rn = pow(LJ_SIG/r,6);
	return -24*LJ_EPS*(2*rn*rn-rn)/(r*r);	
	//return 1.0;	
}

#define loop(x,n) for(int x = 0; x < n; ++x)
#define loop2(x,a,b) for(int x = a; x < b; ++x)

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

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

struct Vector4 {
  double comp[4];
  
  Vector4(double t, double x, double y, double z)
  {
    comp[0] = t; comp[1] = x;  comp[2] = y;  comp[3] = z;
  }
  
  Vector4()
  {
    comp[0] = 0;  comp[1] = 0;  comp[2] = 0;  comp[3] = 0;
  }
};

struct iVector5 {
  int comp[5];
  
  iVector5(int a, int b, int c, int d, int e)
  {
    comp[0] = a;  comp[1] = b;  comp[2] = c;  comp[3] = d;  comp[4] = e;
  }
  
  iVector5()
  {
    comp[0] = 0;  comp[1] = 0;  comp[2] = 0;  comp[3] = 0;  comp[4] = 0;
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

double cross_mag(Vector3 u, Vector3 v)
{
  double w0 = u.comp[1]*v.comp[2] - u.comp[2]*v.comp[1];
  double w1 = u.comp[2]*v.comp[0] - u.comp[0]*v.comp[2];
  double w2 = u.comp[0]*v.comp[1] - u.comp[1]*v.comp[0];
  
  return sqrt(w0*w0 + w1*w1 + w2*w2);
}

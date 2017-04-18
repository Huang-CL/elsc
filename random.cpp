#include"head.h"
extern Ran *rdm;

Ran::Ran(Ullong j): v(4101842887655102017LL), w(1) 
//Constructor. Call with any integer seed (except value of v above).
{
  u = j ^ v; int64();
  v = u; int64();
  w = v; int64();
}

inline Ullong Ran::int64()
{
  u = u * 2862933555777941757LL + 7046029254386353087LL;
  v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
  w = 4294957665U*(w & 0xffffffff) + (w >> 32);
  Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
  return (x + v) ^ w;
}

double Ran::doub()
//Return random double-precision floating value in the range 0. to 1.
{ 
  return 5.42101086242752217E-20 * int64(); 
}

inline Uint Ran::int32()							
//Return 32-bit random integer.
{ return (Uint)int64(); }

double Normal()									
//Return a normal deviate.
{		
  Doub u,v,x,y,q;
  do {
    u = rdm->doub();
    v = 1.7156*(rdm->doub()-0.5);
    x = u - 0.449871;
    y = abs(v) + 0.386595;
    q = pow(x,2) + y*(0.19600*y-0.25472*x);
  } while (q > 0.27597 && (q > 0.27846 || pow(v,2) > -4.*log(u)*pow(u,2)));
  //If u=0,v=0,0.27597<q<0.27846, and log(u) is not a number. Luckly, if one component of A>B expression is not a number, it return false. So this condition will not throw a error.
  return v/u;
}

double Latitude()
{
  double x=1-2*rdm->doub();  //In case of an accidental result
  if(x>1)
    x=1-1e-300;
  else if(x<-1)
    x=-1+1e-300;
  return acos(x);
}

double Longtitude()	//Return x/(2pi)
{
  return 2*pi*rdm->doub();
}

vecdoub::vecdoub(int n)
{
  nn=n;
  x=new double[nn];
}

vecdoub::vecdoub(vecdoub &a)
{
  nn=a.size();
  x=new double[nn];
  for(int i=0;i<nn;i++)
    x[i]=a.x[i];
}

vecdoub::~vecdoub()
{
  if(x!=NULL)
    delete[] x;
  x=NULL;
}

vecdoub::vecdoub(int n, double* a)
{
  nn=n;
  x=new double[nn];
  for(int i=0;i<n;i++)
    x[i]=a[i];
}

void vecdoub::resize(int n)
{
  nn=n;
  delete[] x;
  x=new double[nn];
}

void vecdoub:: operator *= (const double &a)
{
  for(int i=0;i<nn;i++)
    this->x[i]*=a;
}

vecdoub& vecdoub::operator = (const vecdoub &a)
{
  if(this != &a)		// protect against invalid self-assignment
  {
    nn=a.nn;
    double *newx=new double[nn];
    for(int i=0;i<nn;i++)
      newx[i]=a.x[i];
    if(x!=NULL)
      delete[] x;
    x=newx;
  }
  return *this;
}

vec3::vec3(vec3 &a)
{
  for(int i=0;i<3;i++)
    x[i]=a.x[i];
}

vec3::vec3(double* a)
{
  for(int i=0;i<3;i++)
    x[i]=a[i];
}

vec3 vec3:: operator * (const double &a)
{
  double temp[3];
  for(int i=0;i<3;i++)
    temp[i]=x[i]*a;
  vec3 result(temp);
  return result;
}

double vec3::operator * (const vec3 &a)
{
  double product=0.;
  for(int i=0;i<3;i++)
    product+=x[i]*a.x[i];
  return product;
}

vec3 vec3::operator + (const vec3 &a)
{
  double temp[3];
  for(int i=0;i<3;i++)
    temp[i]=x[i]+a.x[i];
  vec3 result(temp);
  return result;
}

vec3 vec3::operator - (const vec3 &a)	
{
  double temp[3];
  for(int i=0;i<3;i++)
    temp[i]=x[i]-a.x[i];
  vec3 result(temp);
  return result;
}

// void photon::save()
// {
//   freq.width(10);
//   freq<<setw(10)<<c*(nu0-nu)/nu0<<'\t';
// }



photon::photon(photon& a)		//Copy constructor function		
{
  x=a.x;
  nnu=a.nnu;
  nu=new double[nnu];
  for(int i=0;i<nnu;i++)
    nu[i]=a.nu[i];
  theta=a.theta;
  phi=a.phi;
}


photon::photon(int nv, double *nu_in, double th, double ph)
{
  x.x[0]=0;x.x[1]=0;x.x[2]=0;
  theta=th; phi=ph;
  nnu=nv;
  nu=new double[nnu];
  for(int i=0;i<nnu;i++)
    nu[i]=nu_in[i];
}

photon& photon::operator = (const photon &a)
{
  if(this != &a)		// protect against invalid self-assignment
  {
    x=a.x;
    theta=a.theta;
    phi=a.phi;
    nnu=a.nnu;
    double *newnu=new double[nnu];
    for(int i=0;i<nnu;i++)
      newnu[i]=a.nu[i];
    if(nu!=NULL)
      delete[] nu;
    nu=newnu;
  }
  return *this;
}

photon::~photon()
{
  if(nu!=NULL)
    delete[] nu;
  nu=NULL;
}

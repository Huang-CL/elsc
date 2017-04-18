#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>	//seed
#include <cstdlib>	//atoi,atof
#include "mpi.h"	// mpi
#include <cstring>	//strcat
#include <sstream>	//stringstream
#include <sys/time.h>		// gettimeofday
#include <vector>

using namespace std;

typedef unsigned long long int Ullong;
typedef double Doub;
typedef unsigned int Uint;

const double pi=3.14159265359;
const double G=6.67384E-8;
const double h=6.62606957E-27;
const double c=2.997924589E10;
const double kb=1.3806488E-16;
const double e=4.80320451E-10;
const double me=9.10938291E-28;
const double mp=1.672621777E-24;
const double sigmaT=0.66524585E-24;
const double Msun=1.9891E33;
const double wmin=1E-10;		//Threshold of ignored photon
const double nu0=1E15;			//Spectral line frequency in
const int Nbin=5000;			// effective bin to save the result (not include overflow bin)


/*Implementation of the highest quality recommended generator. 
  The constructor is called with an integer seed and creates an instance of the generator. 
  The member functions int64, doub, and int32 return the next values in the random sequence, 
  as a variable type indicated by their names. The period of the generator is 3.138*10^57.*/

template<class T>
inline const T &MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

struct overflow{};
//struct used for exception handling. If the array is not large enough
struct onepoint{};
//struct used for exception handling. When interpolate an one element array.

struct vecdoub			//n-D vector
{
  double *x;
  int nn;
  vecdoub(int n=0);
  vecdoub(vecdoub &a);
  vecdoub(int n, double* a);
  void resize(int n);
  int size() const{return nn;}		
  void operator *= (const double &a);			//Vector multiple a number
  vecdoub& operator = (const vecdoub &a);
  ~vecdoub();
};
typedef vecdoub* vecptr;

struct vec3			//3-D vector
{
  vec3(){};
  vec3(vec3 &);
  vec3(double *a);
  vec3 operator * (const double &a);	//Vector multiple a number
  double operator * (const vec3 &a);	//Dot product
  vec3 operator + (const vec3 &a);
  vec3 operator - (const vec3 &a);
  double x[3];
};

struct Ran 
{
  Ullong u,v,w;
  Ran(Ullong j); 
  inline Ullong int64();		
  //Return 64-bit random integer.
  //Inline function. Compiler will copy this part of code to all places that this function was called.
  double doub();								
  //Return random double-precision floating value in the range 0. to 1.
  inline Uint int32();
  //Return 32-bit random integer.
};

double Normal();		//generate normal deviates.

double Latitude();
//generate the theta component of a uniform solid angle deviate

double Longtitude();
//generate the phi component of a uniform solid angle deviate

struct photon
{
  photon(photon& a);		//Copy constructor function		
  photon(int nv=0){nnu=nv;nu=new double[nnu];x.x[0]=0;x.x[1]=0;x.x[2]=0;}
  photon(int nv, double *nu_in, double th, double ph);
  photon& operator = (const photon &a);
  ~photon();
  
  double theta, phi;   //Frequency, transmit direction, escape probablity. 
  int nnu;	       // a photon has nnu frequencies corresponds to nv velocity profiles.
  double *nu;
  //Using w<0 to represent the destroyed photon, w>1 to represent the escaped photon
  vec3 x;		//The coordinate of the photon
};


struct background
{
  int n, nnu;
  double rmin, rmax, mass;
  vecdoub *r, *rho, *emiss, *tau, *emisscum;	
  //Arrays of radius, density, emission rate(erg/vol/time), optical depth, cumulate emission rate
  vecptr *vr, *T;
  //an array of radial velocity array, temperature
  background(double tauT, double power, double abs, double r1, double r2, int nv, double *t, double *vrmin, double *vrmax, double mue=1.18);
  //Constructor function. Need critical optical depth, the power of density profile, absorption probability, destroy radius, escaping radius, number of velocity pairs, array of temperature, minimum & maximum radial velocity and mass of gas per number electron, 1.18 assumes 30% He mass fraction and H, He are fully ionized.
  ~background();
};

struct spectra
{
  int nv;
  double range, binc;  //how many bins that one unit of frequency goes into
  double **bbin, **nbin, **bsumw2, **nsumw2;

  spectra(int nv, double range);
  ~spectra();
  void reset();
  void reduce_sum(MPI::Intracomm worker_comm, int rank);
  /* sum the bins in each processor inside worker communicator and give it to the rank processor */
  void fillnu(double nu, int iv, int ncoll, double increment, double deviate); // fill the frequency into histogram
  void print(string filename, Ullong N, double tauT, double power, double abs, double r1, double r2, double range, int inv, double T, double vrmin, double vrmax, double mass);
};
  

		  
double radius(const vec3 x);

double interp(vecdoub &xv, const vecdoub *yv, double x);	
//Given an array of abscissas xj and function values yi=f(xi), calculate the linear interpolation at x

// double interp(vecdoub &xv, const vecdoub &yv, double x);	
// //Given an array of abscissas xj and function values yi=f(xi), calculate the linear interpolation at x

vec3 interp(vecdoub &xv, const vec3 *yv, double x);
//Simultaneously interpolate multi-vectors.

photon redistribute(double *T, const vec3 *vbar, photon &in);	
//given the temperature of electron gas and input photon data, generate a output photon

photon move2(background *bg, photon &in, double abs, double &escprob, bool &escaped, double mue=1.18);
//Using particle escape weighting method to solve the Monte Carlo simulation. 

//void fitredist(double T,photon in);						//verify the result of redistribute


void mpiwrap(MPI::Intracomm worker_comm,  string *filenames, Ullong N, double tauT, double power, double abs, double r1, double r2, double range, int nv, double *t, double *vrmin, double *vrmax, double mue=1.18);
//total photon number, critical optical depth, the power of density profile, absorption possibility, destroy radius, escaping radius, velocity of last bin, temperature, minimum & maximum radial velocity, mass of gas per number electron, 1.18 assumes 30% He mass fraction and H, He are fully ionized.


void sphere(spectra *spc, background *bg, double abs, double *t, int nv, double mue); // calculate one photon, given a background profile, absorption fraction per scattering, temperatures, number of variables pairs

#include "head.h"
extern ofstream *info;
extern string *bkupfile;
extern Ran *rdm;

background::background(double tauT, double power, double abs, double r1, double r2, int nv, double *t, double *vrmin, double *vrmax, double mue)
{
  double rhoc, fact, norm, kappa=sigmaT/mue/mp/(1-abs);
  nnu=nv;
  n=1000;			//assign space
  rmin=r1;
  rmax=r2;
  vr=new vecptr[nv];
  T=new vecptr[nv];
  r=new vecdoub(n);
  rho=new vecdoub(n);
  emiss=new vecdoub(n);
  tau=new vecdoub(n);
  emisscum=new vecdoub(n);
  for(int i=0;i<nv;i++)
  {
    vr[i]=new vecdoub(n);
    T[i]=new vecdoub(n);
  }
  
  for(int i=0;i<n;i++){			//Initialize arrays
    emiss->x[i]=0;
    r->x[i]=rmin*pow(rmax/rmin,(double)i/(n-1));
    rho->x[i]=pow(rmin/r->x[i],power);
    //rho->x[i]=pow(rmin/r->x[i],2)/(1+pow(r->x[i]/1E15,10)); // The density profile of Model B in Chugai 2001.  Assume the photosphere of the SN is at the Shock front and the photosphere is an absorption bottom boundary.  R_p=R_s=rmin=5E14 cm.  rho \propto (rmin/r)^2/(1+(r/r_tran)^10),  The transition radius is 10^15 cm for model B.  vmin=40 km/s, vmax=1000 km/s.
    emiss->x[i]=/*2.3E21*/rho->x[i]*rho->x[i];
    for(int j=0;j<nv;j++)
    {
      T[j]->x[i]=t[j];
      vr[j]->x[i]=vrmin[j]*1E5+vrmax[j]*1E5*pow(rmin/r->x[i],2.0);
    }
  }
  tau->x[n-1]=0;		//Integrate the optical depth
  for(int i=n-2;i>=0;i--)
    tau->x[i]=tau->x[i+1]+(r->x[i+1]-r->x[i])*rho->x[i+1]*kappa;
  rhoc=tauT/tau->x[0];		//Calculate the normalization factor
  *rho*=rhoc;	//Normalize the density, optical depth and emission rate
  *tau*=rhoc;
  *emiss*=(rhoc*rhoc);
  emisscum->x[0]=0;	//Integrate the cumulated emission and normalize to 1.
  for(int i=1;i<n;i++)
    emisscum->x[i]=emisscum->x[i-1]+/*16*pow(pi,2)/3**/(pow(r->x[i],3.0)-pow(r->x[i-1],3.0))*emiss->x[i-1];
  norm=emisscum->x[n-1];
  *emisscum*=(1./norm);
  mass=0;			//Integrate the mass
  fact=4*pi/3;
  for(int i=0;i<n-1;i++)
    mass+=fact*(pow(r->x[i+1],3.0)-pow(r->x[i],3.0))*rho->x[i];
}

background::~background()
{

  delete r;
  delete rho;
  delete emiss;
  delete tau;
  delete emisscum;
  r=rho=emiss=tau=emisscum=NULL;
  for(int i=0;i<nnu;i++)
  {
    delete vr[i];
    delete T[i];
  }
  delete[] vr;
  delete[] T;
}


void sphere(spectra *spc, background *bg, double abs, double *t, int nv, double mue) // calculate one photon, given a background profile, absorption fraction per scattering, temperatures, number of variables pairs.
{
  vec3* vbar=new vec3[nv];           //velocity of medium
  photon data(nv);
  double theta, phi, rr, r, vr;		//Variables to produce the location of photon
  double w, escprob;	//The weight of photon
  double increment, deviate;
  double delta;		//Factor of frequency change by change reference frame
  bool escaped;
  vec3 v;		//Peculiar velocity
  vec3 l;		//Photon direction in Cartesian coordinates
  double *vtherm=new double[nv];//proton velocity mean square root in one direction
  for(int i=0;i<nv;i++)
    vtherm[i]=sqrt(kb*t[i]/mp);

  // rr=((double)(n+1))/((double)(N+2));
  rr=rdm->doub();
  r=interp(*bg->emisscum,bg->r,rr);
  int ncoll=0;
  phi=Longtitude();		//Randomly generate a photon
  theta=Latitude();

  data.x.x[0]=r*sin(theta)*cos(phi);
  data.x.x[1]=r*sin(theta)*sin(phi);
  data.x.x[2]=r*cos(theta);
  data.phi=Longtitude();	//The moving direction of photon
  data.theta=Latitude();
  l.x[0]=sin(data.theta)*cos(data.phi);
  l.x[1]=sin(data.theta)*sin(data.phi);
  l.x[2]=cos(data.theta);

  for(int i=0;i<nv;i++)
  {
    vr=interp(*bg->r,bg->vr[i],r);
    vbar[i]=data.x*(vr/r);//The velocity of medium in Cartesian coordinates
    v.x[0]=vbar[i].x[0]+vtherm[i]*Normal();
    //generate a random velocity in one direction
    v.x[1]=vbar[i].x[1]+vtherm[i]*Normal();
    v.x[2]=vbar[i].x[2]+vtherm[i]*Normal();
    delta=l*v;
    data.nu[i]=nu0*(1+delta/c);
  }

//nu0 is the frequency in medium comoving frame (lab frame), nu is in observer's frame
  w=1.0;
  do
  {
    data=move2(bg,data,abs,escprob,escaped,mue);
    if(escaped)
    {
      //don't have to save the photon that contribute very little to the overall shape
      increment=w*escprob;
      deviate=increment*increment;
      for(int i=0;i<nv;i++)
	spc->fillnu(data.nu[i],i,ncoll,increment,deviate);
    }
    ncoll++;
    w=w*(1-escprob)*(1-abs);
    //For each scattering, there are abs of the probability that the photon absorbed by the medium. w is the weight that this photon haven't escaped or absorbed.
    r=radius(data.x);
    for(int i=0;i<nv;i++)
    {
      vr=interp(*bg->r,bg->vr[i],r);
      vbar[i]=data.x*(vr/r);//The velocity of medium in Cartesian coordinates
    }
    data=redistribute(t,vbar,data);
  }while(w>1E-7);
  delete[] vbar;
  delete[] vtherm;
  vbar=NULL;
}
  

void mpiwrap(MPI::Intracomm worker_comm, string *filenames, Ullong N, double tauT, double power, double abs, double r1, double r2, double range, int nv, double *t, double *vrmin, double *vrmax, double mue)
//total photon number, critical optical depth, the power of density profile, absorption possibility, destroy radius, escaping radius, binning range, temperature, minimum & maximum radial velocity, mass of gas per number electron
{
  int size = MPI::COMM_WORLD.Get_size();
  int myid = MPI::COMM_WORLD.Get_rank();
  int tag_request=1,tag_command=2;
  int nworkers = size-1;
  int nshutdown=0;
  int masterServer=size-1;
  bool request=true;
  MPI::Status status;
  Ullong pcount=myid;

  spectra *spc;
  background *bg;
  
  if(myid==masterServer)
  {
    pcount=nworkers;
    if(pcount>N)		// more workers than the number of photon
      nshutdown=nworkers-N;

    while(pcount<N || nshutdown<nworkers)
    {
      MPI::COMM_WORLD.Recv(&request,1,MPI::BOOL,MPI::ANY_SOURCE,tag_request,status);
      if(pcount<N)
      {
	MPI::COMM_WORLD.Send(&pcount,1,MPI::UNSIGNED_LONG_LONG,status.Get_source(),tag_command);
	//my progress bar
	if(100*(pcount-nworkers+1)/N>100*(pcount-nworkers)/N)
	  cout<<100*(pcount-nworkers+1)/N<<'%'<<endl;
	pcount++;
      }
      else
      {
	MPI::COMM_WORLD.Send(&pcount,1,MPI::UNSIGNED_LONG_LONG,status.Get_source(),tag_command);
	//my progress bar
	if(nshutdown==nworkers-1) // in case the percentage is larger than 100%
	  cout<<"100%"<<endl;
	else if(100*(N-nworkers+nshutdown+1)/N>100*(N-nworkers+nshutdown)/N)
	  cout<<100*(N-nworkers+nshutdown+1)/N<<'%'<<endl;
	nshutdown++;
      }
    }
    //cout<<"masterServer finished it's job"<<endl;
  }

  // worker processor
  else
  {
    spc=new spectra(nv,range);
    bg=new background(tauT,power,abs,r1,r2,nv,t,vrmin,vrmax,mue);
    
    while(pcount<N)
    {
      sphere(spc, bg, abs, t, nv, mue);
      MPI::COMM_WORLD.Send(&request,1,MPI::BOOL,masterServer,tag_request);
      MPI::COMM_WORLD.Recv(&pcount,1,MPI::UNSIGNED_LONG_LONG,masterServer,tag_command);
    }
    spc->reduce_sum(worker_comm,0);

    if(myid==0)
      for(int j=0;j<nv;j++)
	spc->print(filenames[j],N,tauT,power,abs,r1,r2,range,j,t[j],vrmin[j],vrmax[j],bg->mass);

    //cout<<"processor "<<myid<<" is ready to clean up"<<endl;  
    // clean up
    delete spc;
    delete bg;
  }
  //cout<<"processor "<<myid<<" finished one mpiwrap"<<endl;  
}

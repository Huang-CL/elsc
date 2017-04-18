#include"head.h"
extern Ran *rdm;

photon redistribute(double *T, const vec3 *vbar, photon &in)	
//given the temperature of electron gas and input photon data, generate a output photon
{
  int nv=in.nnu;
  double *vtherm=new double[nv];	//electron velocity mean square root in one direction

  vec3 v;
  photon out=in;	//the output photon data
  double cosTheta,p,ran, delta;
  do{
    out.theta=Latitude();   //get the output direction
    out.phi=Longtitude();
    cosTheta=sin(in.theta)*sin(out.theta)*cos(in.phi-out.phi)+cos(in.theta)*cos(out.theta);	//calculate the angle between input photon and output photon
    p=0.5*(1+pow(cosTheta,2));
    ran=rdm->doub();
  }while(ran>=p);	//use rejection method to pick out good data point

  for(int i=0;i<nv;i++)
  {
    vtherm[i]=sqrt(kb*T[i]/me);  	//electron velocity mean square root in one direction
    v.x[0]=vbar[i].x[0]+vtherm[i]*Normal();//generate a random velocity in one direction
    v.x[1]=vbar[i].x[1]+vtherm[i]*Normal();
    v.x[2]=vbar[i].x[2]+vtherm[i]*Normal();

    delta=v.x[0]*(sin(out.theta)*cos(out.phi)-sin(in.theta)*cos(in.phi))+v.x[1]*(sin(out.theta)*sin(out.phi)-sin(in.theta)*sin(in.phi))+v.x[2]*(cos(out.theta)-cos(in.theta));
    out.nu[i]=in.nu[i]*(1+delta/c);	
    //calculate the frequency of output photon in choosen direction.
    //cout<<vtherm<<"\t "<<v.x[0]<<"\t "<<v.x[1]<<"\t "<<v.x[2]<<"\t "<<in.x.x[0]<<"\t "<<in.x.x[1]<<"\t "<<in.x.x[2]<<"\t "<<in.nu<<"\t "<<out.nu<<"\t "<<in.theta<<"\t "<<in.phi<<"\t "<<out.theta<<"\t "<<out.phi<<endl;
  }
  delete [] vtherm;
  return out;
}

double radius(const vec3 x){	    //get radial distance
  return sqrt(x.x[0]*x.x[0]+x.x[1]*x.x[1]+x.x[2]*x.x[2]);}

double density(vec3 x, background *bg)
{
  double r=radius(x);
  if(r<0.9999*bg->rmin)	//keep working only when the photon stay in the chosen area.
    cout<<"destroyed"<<endl;
  else if(r>1.0001*bg->rmax)
    cout<<"escaped"<<endl;
  return interp(*bg->r,bg->rho,r);	//get gas density at this point
}

void odeint(background *bg, photon &in, double abs, vec3 n, vec3 *xp, vecdoub &taup, double mue)	
//The function to do the ordinate, build up the array of location and optical depth along the trajectory of photon.
{
  double kappa=sigmaT/mue/mp/(1-abs);
  double step, h;	// step length, 1/3 of the step length
  double rho0, rho1, rho2, rho3;	//density at two ends of the range and two 1/3 mid points.
  double cosTheta,r0;		//the angle between photon trajectory and radial direction
  int i;
  vec3 dx;
  xp[0]=in.x;		//Initialize the first point
  taup.x[0]=0;
  vec3 x0;   			//The coordinate of a point
  double r=radius(xp[0]);
  rho0=density(xp[0],bg);	//initialize the density and radius.
  h=0.001*r;	//first step size. Try to avoid if condition in the loop
  step=h*3.;
  dx=n*h;
  xp[1]=xp[0]+dx*3.;
  r=radius(xp[1]);
  for(i=1;r<bg->rmax&&r>bg->rmin;i++){//Get the optical depth at every point

    x0=xp[i-1]+dx;		//1/3 point. 
    rho1=density(x0,bg);

    x0=x0+dx;			//2/3 point
    rho2=density(x0,bg);

    rho3=density(xp[i],bg);	//right end of the range

    taup.x[i]=taup.x[i-1]+step*(rho0+3*(rho1+rho2)+rho3)*kappa/8;

    rho0=rho3;	//the density of the left end of the range for next loop

    h=0.001*r;
    //Step size is changing. Smaller step is used at small r in order to increase accuracy. 
    step=h*3.;
    dx=n*h;
    xp[i+1]=xp[i]+dx*3.;			//Next point location.
    r=radius(xp[i+1]);
    //The radius of next point. Used to decide next step length and used as a condition to continue the loop.
 
  }

  r0=radius(xp[i-1]);
  cosTheta=n*xp[i-1]/r0;

  if(r>=bg->rmax&&cosTheta<0.1)
    //The last step cannot go out of the gas region. If cos(theta) is too small, the approximation doesn't work.
    h=0.33*bg->rmax*(sqrt(cosTheta*cosTheta+2*(bg->rmax-r0)/r0)-cosTheta);
  else if(r>=bg->rmax)
    h=0.33*(bg->rmax-r0)/cosTheta;
  else if(r<=bg->rmin&&cosTheta>-0.1)
    h=-0.33*bg->rmin*(cosTheta+sqrt(cosTheta*cosTheta-2*(r0-bg->rmin)/r0));
  else if(r<=bg->rmin)
    h=0.33*(bg->rmin-r0)/cosTheta;

  dx=n*h;
  xp[i]=xp[i-1]+dx*3.;
  step=h*3;

  x0=xp[i-1]+dx;		//1/3 point. 
  rho1=density(x0,bg);

  x0=x0+dx;			//2/3 point
  rho2=density(x0,bg);

  rho3=density(xp[i],bg);	//right end of the range

  taup.x[i]=taup.x[i-1]+step*(rho0+3*(rho1+rho2)+rho3)*kappa/8;
  i++;
  taup.nn=i;
  if(i>10000)
    throw overflow();
}

double tauscatfunc(double r, double tau)					//return tauscat=-ln[1-r*(1-exp(-tau0))]
{
  if(tau>1E-5)			//See ref 1983 P313
    return -log(1.0-r*(1.0-exp(-tau)));
  else		
    //When tau is very small, the cancellation error in 1-exp(-tau) will largely decrease the accuracy of result
    return r*tau+(pow(r,2)-r)*pow(tau,2)/2+(r-3*pow(r,2)+2*pow(r,3))*pow(tau,3)/6;
  //Using Taylor expension to the third order of tau to solve it.
}

photon move2(background *bg, photon &in, double abs, double &escprob, bool &escaped, double mue)
{
  int kmax=10000;
  //the maximun number of points used to calculate the optical depth
  photon out=in;
  //cout<<out.nu[0]<<'\t'<<out.nnu<<'\t'<<out.theta<<'\t'<<in.nu[0]<<'\t'<<in.nnu<<'\t'<<in.theta<<endl;
  vec3 n;	//transmit direction vector
  n.x[0]=sin(in.theta)*cos(in.phi);n.x[1]=sin(in.theta)*sin(in.phi);n.x[2]=cos(in.theta);
  double AA=in.x*n;
  //The distance from current location to the point that has smallest r alone the transit path
  //x(l)=x0+n l; Both sides dot itself. r^2(l)=r0^2+2 x0*n l+l^2=r0sq+2 AA l+l^2
  double r0sq=in.x*in.x;
  double bsq=r0sq-pow(AA,2);							//squre of impact parameter
  if(AA>=0||bsq>pow(bg->rmin,2))
    escaped=true;								
  //If escaped=true, this photon may escaped if it is not scattered (impossible to be distroyed). If escpaed=false, this photon may distroyed
  else
    escaped=false;
  vec3 *xp=new vec3[kmax];	// a coordinate array along the line
  vecdoub *taup=new vecdoub(kmax);						
  //The arrays to store the data along the trajectory
  odeint(bg, in, abs, n, xp, *taup,mue);						
  //Initialize the optical depth array, used for interpolation.
  int ng=taup->size();
  double taubase=taup->x[ng-1];		       					
  //The total optical depth along this trajectory
  if(taubase<1E-4){
    //if the input already very close to the boundary, and will move out the medium, just set it certainly will escaped or destroyed
    escprob=1;
    delete taup;
    delete[] xp;
    return out;
  }
  double r1=rdm->doub();							
  //Produce a random number to indicate the position that photon is scattered if it neither escape nor destory.
  double tauscat=tauscatfunc(r1,taubase), rout;					
  //The photon is scattered after it travel the medium which optical depth is tauscat
  out.x=interp(*taup,xp,tauscat);	//Find the location to be scattered
  rout=radius(out.x);
  //cout<<radius(in.x)<<"\t "<<taubase<<"\t "<<r1<<"\t "<<tauscat<<"\t "<<rout<<endl;
  if(rout<bg->rmin&&!escaped)
    //If the photon move out of the region, put it back to the medium
    out.x=(out.x*(1/rout))*1.00001*bg->rmin;
  else if(rout>bg->rmax&&escaped)
    out.x=(out.x*(1/rout))*0.99999*bg->rmax;

  escprob=exp(-min(taubase,100.0));
  //The probability that the photon still remain in the medium region. To avoid some extreme cases. If taubase>100, we definitely will ignore this photon
  delete taup;
  delete[] xp;
  return out;
}




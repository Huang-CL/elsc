#include "head.h"

spectra::spectra(int nv_in, double range_in):nv(nv_in),range(range_in)
{
  binc=Nbin*c*5E-6/nu0/range;
  bbin=new double*[nv];
  nbin=new double*[nv];
  bsumw2=new double*[nv];
  nsumw2=new double*[nv];	//binning the narrow component and broad component seperately, sumw2 stores the summation of weight square used to calculate the error. The error per bin will be computed as sqrt(sum of squares of weight) for each bin.

  for(int j=0;j<nv;j++)
  {
    bbin[j]=new double[Nbin+2];
    nbin[j]=new double[Nbin+2];
    bsumw2[j]=new double[Nbin+2];
    nsumw2[j]=new double[Nbin+2];
    for(int i=0;i<Nbin+2;i++)
      bbin[j][i]=nbin[j][i]=bsumw2[j][i]=nsumw2[j][i]=0;
  }
}

spectra::~spectra()
{
  for(int j=0;j<nv;j++)
  {
    delete [] bbin[j];
    delete [] nbin[j];
    delete [] bsumw2[j];
    delete [] nsumw2[j];
  }
  delete [] bbin;
  delete [] nbin;
  delete [] bsumw2;
  delete [] nsumw2;
}

void spectra::reset()
{
  for(int j=0;j<nv;j++)
    for(int i=0;i<Nbin+2;i++)
      bbin[j][i]=nbin[j][i]=bsumw2[j][i]=nsumw2[j][i]=0;
}

void spectra::fillnu(double nu, int iv, int ncoll, double increment, double deviate) // fill the frequency into histogram
{
  int k=int((nu0-nu)*binc+Nbin/2+1);
  if(k<0)
    k=0;
  else if(k>Nbin+1)
    k=Nbin+1;
  if(!ncoll)
  {
    nbin[iv][k]+=increment;
    nsumw2[iv][k]+=deviate;
  }
  else
  {
    bbin[iv][k]+=increment;
    bsumw2[iv][k]+=deviate;
  }
}

void spectra::reduce_sum(MPI::Intracomm worker_comm, int rank)
/* sum the bins in each processor inside worker communicator and give it to the rank processor */
{
  int myid = MPI::COMM_WORLD.Get_rank();

  if(myid==rank)
  {
    for(int i=0;i<1;i++)
    {
      worker_comm.Reduce(MPI::IN_PLACE, bbin[i],Nbin+2,MPI::DOUBLE,MPI::SUM,rank);
      worker_comm.Reduce(MPI::IN_PLACE, nbin[i],Nbin+2,MPI::DOUBLE,MPI::SUM,rank);
      worker_comm.Reduce(MPI::IN_PLACE, bsumw2[i],Nbin+2,MPI::DOUBLE,MPI::SUM,rank);
      worker_comm.Reduce(MPI::IN_PLACE, nsumw2[i],Nbin+2,MPI::DOUBLE,MPI::SUM,rank);
    }
  }
  else
  {
    for(int i=0;i<1;i++)
    {
      worker_comm.Reduce(bbin[i], bbin[i],Nbin+2,MPI::DOUBLE,MPI::SUM,rank);
      worker_comm.Reduce(nbin[i], nbin[i],Nbin+2,MPI::DOUBLE,MPI::SUM,rank);
      worker_comm.Reduce(bsumw2[i], bsumw2[i],Nbin+2,MPI::DOUBLE,MPI::SUM,rank);
      worker_comm.Reduce(nsumw2[i], nsumw2[i],Nbin+2,MPI::DOUBLE,MPI::SUM,rank);
    }
  }
}

void spectra::print(string filename, Ullong N, double tauT, double power, double abs, double r1, double r2, double range, int inv, double T, double vrmin, double vrmax, double mass)
{
  double narrow=0, total=0;
  ofstream info(filename.c_str(),ios_base::app);
  if(!info)
  {
    cout<<"fail to open result file "<<filename<<endl;
  }
  
  info<<"r^-"<<power<<" density profile \t absorption="<<abs*100<<"% \t # of photon="<<N<<"\t tau="<<tauT<<"\t shell radius:"<<r1<<'/'<<r2<<"cm\t velocity, v=vmin+vmax*(r1/r)^2:"<<vrmin<<'/'<<vrmax<<"km/s \t total mass:"<<mass<<" g"<<endl;
  info<<T<<" K temperature"<<endl; 
  info<<range<<"\t range in km/s"<<endl;
  for(int i=0;i<Nbin+2;i++)
  {
    info<<i<<'\t'<<bbin[inv][i]<<'\t'<<bsumw2[inv][i]<<'\t'<<nbin[inv][i]<<'\t'<<nsumw2[inv][i]<<endl;
    total+=bbin[inv][i]+nbin[inv][i];
    narrow+=nbin[inv][i];
  }
  info<<"narrow \t total"<<endl<<narrow<<"\t"<<total<<endl;

  info.close();
}


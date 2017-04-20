#include "head.h"
string *filenames;		//the string to save the result of the data

struct timeval start_time;
Ran *rdm;

int main(int argc, char *argv[])
{
  if(argc!=2)
    argv[1]=(char *)"arguments.txt";
  ifstream arg(argv[1]);
  if(!arg){
    cout<<"wrong file"<<endl;
    return 0;
  }
  string line,buf;
  vector<string> value;
  int counts;
  string datafile,headfile;
  stringstream ss;
  double *t, *vrmin, *vrmax;
  int nv;

  int masterServer;
  int ranks[1];
  MPI::Group world_group;
  MPI::Group worker_group;
  MPI::Intracomm worker_comm;
  MPI::Init();
  int size = MPI::COMM_WORLD.Get_size();
  //int myid = MPI::COMM_WORLD.Get_rank();
  
  //Build a worker communicator in order to use MPI::Reduce
  world_group = MPI::COMM_WORLD.Get_group();
  masterServer=size-1;		// better not change this setting. Or change setting in mpi as well.
  ranks[0]=masterServer;
  worker_group=world_group.Excl(1,ranks); // exclude master processor from the communicator
  worker_comm=MPI::COMM_WORLD.Create(worker_group);
  worker_group.Free();
  world_group.Free();


  gettimeofday(&start_time,NULL);
   // microsecond has 1 000 000
  // Assuming you did not need quite that accuracy
  // Also do not assume the system clock has that accuracy.
  rdm = new Ran(start_time.tv_sec * 1000000 + start_time.tv_usec);
  // If you use 100 (rather than 1000) the seed repeats every 248 days.

  // Do not make the MISTAKE of using just the tv_usec
  // This will mean your seed repeats every second.
  while(getline(arg,line))
  {
    counts=0;
    ss.str(std::string());
    ss.clear();
    ss<<line;
    while(ss >> buf)
    {
      value.push_back(buf);
      counts++;
    }
    if (counts != 9)
    {
      cout<<"counts="<<counts<<"\n file path, N, tau, power, abs, r1, r2, range, nv"<<endl;
      value.clear();
      getline(arg,line);
      continue;
    }
    else
    {
      getline(arg,line);
      nv=atoi(value[8].c_str());
      counts=0;
      ss.str(std::string());
      ss.clear();
      ss<<line;
      while(ss >> buf)
      {
    	value.push_back(buf);
    	counts++;
      }
      if(counts!=3*nv)
      {
    	cout<<"counts="<<counts<<"\t we need "<<nv<<" pairs of temperature, min velocity and max velocity"<<endl;
    	value.clear();
    	continue;
      }
      else
      {
      	t=new double[nv];
      	vrmin=new double[nv];
      	vrmax=new double[nv];
	filenames=new string[nv];
      	for(int i=0;i<nv;i++)
      	{
      	  filenames[i]=value[0]+value[9+i*3]+"_"+value[10+i*3]+"_"+value[11+i*3]+".info";
      	  t[i]=atof(value[9+i*3].c_str());
      	  vrmin[i]=atof(value[10+i*3].c_str());
      	  vrmax[i]=atof(value[11+i*3].c_str());
      	}
      	mpiwrap(worker_comm, filenames, atof(value[1].c_str()),atof(value[2].c_str()),atof(value[3].c_str()),atof(value[4].c_str()),atof(value[5].c_str()),atof(value[6].c_str()),atof(value[7].c_str()),nv,t,vrmin,vrmax);

      	value.clear();
      	delete[] t;
      	delete[] vrmin;
      	delete[] vrmax;
	delete[] filenames;
      }
    }
    //cout<<"processor "<<myid<<" finished one calculation"<<endl;  
  }

  //cout<<"processor "<<myid<<" is ready to finalize"<<endl;  
  MPI::Finalize();
  
  return 0;
}

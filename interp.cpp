#include"head.h"

double interp(vecdoub &xv, const vecdoub *yv, double x)
{
  int n=xv.size();
  int ju,jm,jl,j;
  if(n<2)throw("locate size error");
  bool ascnd=(xv.x[n-1] >= xv.x[0]);						
  //True if ascending order of table, false otherwise.
  jl=0;		//Initialize lower
  ju=n-1;	//and upper limits.
  while (ju-jl > 1) {//If we are not yet done,
    jm = (ju+jl) >> 1;//compute a midpoint, (ju+jl)/2
    if ((x >= xv.x[jm]) == ascnd)
      jl=jm;	//and replace either the lower limit
    else
      ju=jm;	//or the upper limit, as appropriate.
  }		//Repeat until the test condition is satisfied.
  j=MAX(0,MIN(n-2,jl));
  //Even if the x is larger than the maximum value or smaller than the minimum value, the code won't crash. And the possible largest result is no larger than n-2
  if(xv.x[j]==xv.x[j+1]) 
    return yv->x[j];	//Table is defective, but we can recover.
  else 
    return yv->x[j] + ((x-xv.x[j])/(xv.x[j+1]-xv.x[j]))*(yv->x[j+1]-yv->x[j]);
}


// double interp(vecdoub &xv, const vecdoub &yv, double x)
// {
//   int n=xv.size();
//   int ju,jm,jl,j;
//   if(n<2)throw("locate size error");
//   bool ascnd=(xv.x[n-1] >= xv.x[0]);						
//   //True if ascending order of table, false otherwise.
//   jl=0;		//Initialize lower
//   ju=n-1;	//and upper limits.
//   while (ju-jl > 1) {//If we are not yet done,
//     jm = (ju+jl) >> 1;//compute a midpoint, (ju+jl)/2
//     if ((x >= xv.x[jm]) == ascnd)
//       jl=jm;	//and replace either the lower limit
//     else
//       ju=jm;	//or the upper limit, as appropriate.
//   }		//Repeat until the test condition is satisfied.
//   j=MAX(0,MIN(n-2,jl));
//   //Even if the x is larger than the maximum value or smaller than the minimum value, the code won't crash. And the possible largest result is no larger than n-2
//   if(xv.x[j]==xv.x[j+1]) 
//     return yv.x[j];	//Table is defective, but we can recover.
//   else 
//     return yv.x[j] + ((x-xv.x[j])/(xv.x[j+1]-xv.x[j]))*(yv.x[j+1]-yv.x[j]);
// }



vec3 interp(vecdoub &xv, const vec3 *yv, double x)
{
  vec3 out;
  int n=xv.size();
  int ju,jm,jl,j;
  if(n<2){
    out=*yv;
    cout<<"n="<<n<<endl;
    //  throw("locate size error");
    return out;
  }
  bool ascnd=(xv.x[n-1] >= xv.x[0]);						
  //True if ascending order of table, false otherwise.
  jl=0;	//Initialize lower
  ju=n-1;	//and upper limits.
  while (ju-jl > 1) {	//If we are not yet done,
    jm = (ju+jl) >> 1;	//compute a midpoint, (ju+jl)/2
    if ((x >= xv.x[jm]) == ascnd)
      jl=jm;		//and replace either the lower limit
    else
      ju=jm;		//or the upper limit, as appropriate.
  }		       	//Repeat until the test condition is satisfied.
  j=MAX(0,MIN(n-2,jl));
  int m;
  if(xv.x[j]==xv.x[j+1])//Table is defective, but we can recover. 
    for(m=0;m<3;m++)
      out.x[m]=yv[m].x[j];
  else 
    for(m=0;m<3;m++)
      out.x[m]=yv[j].x[m] + ((x-xv.x[j])/(xv.x[j+1]-xv.x[j]))*(yv[j+1].x[m]-yv[j].x[m]);
  return out;
}


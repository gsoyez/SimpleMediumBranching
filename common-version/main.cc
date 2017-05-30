#include <iostream>
#include "event.hh"
#include "generator.hh"
#include <iomanip>
#include <gsl/gsl_math.h>
using namespace std;


int main(){

  cout << "Running"<<endl;

  double testTime=0.1;
  double epsilon=1e-3;
  GeneratorInMedium gen;
  gen.generateEvent(testTime,epsilon);

  const vector<Particle> &testvec=gen.event().particles();

  cout <<setw(8) << "index" <<setw(8)<< "parent"<<setw(8) <<"child1"<<setw(8) <<"child2";
  cout <<setw(8) <<"start"<<setw(8) <<"end"<<setw(8) <<"x"<<setw(8)<<"final"<<endl;
  cout << setprecision(2);

  for (unsigned int i = 0; i < testvec.size(); i++ ) {
    cout <<setw(8) << i<<setw(8)<< testvec[i].parent();
    cout <<setw(8) << testvec[i].child1();
    cout <<setw(8) << testvec[i].child2();
    cout <<setw(8) << testvec[i].startTime();
    cout <<setw(8) << testvec[i].endTime();
    cout <<setw(8) << testvec[i].x();
    cout <<setw(8) << testvec[i].is_final();
    cout <<endl;
  }

  cout << "Final x-values:"<<endl;
  const vector<double> final_x=gen.event().final_particles();
  for (unsigned int i = 0; i < final_x.size(); i++ ) {
    cout << final_x[i] <<endl;
  }







  ///Checking bin close to 0.5 and first bin in time
  cout << "Checking fluctuations, takes less than a minute at current number of iterations" <<endl;
  GeneratorInMedium v;
  double cutoff = 1e-10;
  v.set_cutoff(cutoff);
  double alpha = (4-8*cutoff)/sqrt(cutoff*(1-cutoff));///<integral of simplified kernel
  unsigned int iter1=1e3;///usually 1e3 is good
  unsigned int iter2=1e4;///<adjust to follow behaviour

  double totalArea=399993.0/2;
  double expectedZ=0.0389833*iter2/totalArea;
  double expectedT=0.000199977*iter2;
  double deviationZ=0;
  double deviationT=0;

  double binZ,binT;
  for (unsigned int i=1; i<iter1;i++){
    binZ=0;
    binT=0;
      for (unsigned int j=1; j<iter2;j++){
        v.generateBranching(1);
        double z=v.z();
        double t=v.t();
          if (z>0.49&&z<0.5){
            binZ+=1;
          }
          if (t<10*cutoff){
            binT+=1;
          }
      }
      deviationZ += (binZ-expectedZ)*(binZ-expectedZ);
      deviationT += (binT-expectedT)*(binT-expectedT);

  }

  cout << setprecision(4) <<"deviationZ : " << sqrt(deviationZ/iter1) << " root n " << sqrt(expectedZ) << endl;
  cout << setprecision(4) <<"deviationT : " << sqrt(deviationT/iter1) << " root n " << sqrt(expectedT) << endl;










  return 0;
}

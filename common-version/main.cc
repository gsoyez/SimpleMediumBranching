#include <iostream>
#include <fstream>
#include <iomanip>
#include "event.hh"
#include "generator.hh"
#include <iomanip>
#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
using namespace std;


int main(){

  cout << "Running"<<endl;

  double test_time=0.1;
  double epsilon=1e-3;

  GeneratorInMedium gen;
  gen.set_seed(1);
  gen.generate_event(test_time,epsilon);
  const vector<Particle>& testvec=gen.event().particles();

  cout <<setw(8) << "index" <<setw(8)<< "parent"<<setw(8) <<"child1"<<setw(8) <<"child2";
  cout <<setw(8) <<"start"<<setw(8) <<"end"<<setw(8) <<"x"<<setw(8)<<"final"<<endl;
  cout << setprecision(2);

  for (unsigned int i = 0; i < testvec.size(); i++ ) {
    cout <<setw(8) << i<<setw(8)<< testvec[i].parent();
    cout <<setw(8) << testvec[i].child1();
    cout <<setw(8) << testvec[i].child2();
    cout <<setw(8) << testvec[i].start_time();
    cout <<setw(8) << testvec[i].end_time();
    cout <<setw(8) << testvec[i].x();
    cout <<setw(8) << testvec[i].is_final();
    cout <<endl;
  }

  cout << "Final x-values:"<<endl;
  const vector<double> final_x=gen.event().final_particles();
  for (unsigned int i = 0; i < final_x.size(); i++ ) {
    cout << final_x[i] <<endl;
  }



//Should test spectrum against known solution...

  GeneratorInMediumSimple gen2;
  gen2.generate_event(test_time,epsilon);
  gen.set_seed(1);
  const vector<Particle>& testvec2=gen2.event().particles();

  cout <<setw(8) << "index" <<setw(8)<< "parent"<<setw(8) <<"child1"<<setw(8) <<"child2";
  cout <<setw(8) <<"start"<<setw(8) <<"end"<<setw(8) <<"x"<<setw(8)<<"final"<<endl;
  cout << setprecision(2);

  for (unsigned int i = 0; i < testvec2.size(); i++ ) {
    cout <<setw(8) << i<<setw(8)<< testvec2[i].parent();
    cout <<setw(8) << testvec2[i].child1();
    cout <<setw(8) << testvec2[i].child2();
    cout <<setw(8) << testvec2[i].start_time();
    cout <<setw(8) << testvec2[i].end_time();
    cout <<setw(8) << testvec2[i].x();
    cout <<setw(8) << testvec2[i].is_final();
    cout <<endl;
  }

  cout << "Final x-values:"<<endl;
  const vector<double> final_x2=gen2.event().final_particles();
  for (unsigned int i = 0; i < final_x2.size(); i++ ) {
    cout << final_x2[i] <<endl;
  }



//Code below here should move to a separate file later




  ///Plotting histograms
  test_time=0.1;
  epsilon=1e-5;

  /*
  gsl_histogram * h_simple = gsl_histogram_alloc (1000);
  gsl_histogram_set_ranges_uniform (h_simple, epsilon, 1-epsilon);
  */

  //Setting logarithmic bins
  unsigned int bins=100;
  double lowest=-log10(epsilon)/bins;
  double range[bins+1]={};
  for(unsigned int i_range=0;i_range<bins+1;++i_range){
    range[bins-i_range]=1/pow(10,lowest*(i_range));
  }


  gsl_histogram * h_simple = gsl_histogram_alloc (bins);
  gsl_histogram_set_ranges (h_simple, range, bins+1);

  gsl_histogram2d * h_simple2d = gsl_histogram2d_alloc (bins,bins);
  gsl_histogram2d_set_ranges (h_simple2d, range, bins+1,range, bins+1);

  unsigned int iterations = 1e3;
  GeneratorInMediumSimple g;
  for (unsigned int i_it=1; i_it<iterations;i_it++){
    g.generate_event(test_time,epsilon);
    const vector<double> finals=g.event().final_particles();
    for (unsigned int i_x = 0; i_x < finals.size(); i_x++ ) {
      gsl_histogram_increment (h_simple, finals[i_x]);
      for (unsigned int i_y = 0; i_y < finals.size(); i_y++ ) {
        if (i_x!=i_y){  //Checking that the gluons are not the same one
          gsl_histogram2d_increment(h_simple2d,finals[i_x],finals[i_y]);
        }
      }
    }
  }

  ofstream ostr("Results/histogram_simple.dat");
  double binwidth, minval, maxval, meanval, value;
  ostr << "# histogram of simple kernel" << endl;
  ostr << "# columns are min, max, value" << endl;
  for (unsigned int i_h=0; i_h<gsl_histogram_bins(h_simple); ++i_h){
    gsl_histogram_get_range(h_simple, i_h, &minval, &maxval);
    binwidth = maxval-minval;
    meanval=0.5*(minval+maxval);
    meanval=(minval);
    value =sqrt(meanval)*meanval*gsl_histogram_get(h_simple, i_h)/(binwidth*iterations);
    ostr << minval << " " << maxval << " " << value << endl;
  }
  ostr.close();



  ofstream ostr2("Results/histogram2d_simple.dat");
  ofstream ostr3("Results/histogram2d_diag_simple.dat");
  double xbinwidth,ybinwidth, xminval, xmaxval,yminval,ymaxval, xmeanval,ymeanval, value2;
  ostr2 << "# histogram of simple kernel" << endl;
  ostr2 << "# columns are X, Y, Z" << endl;
  ostr << "# histogram of simple kernel" << endl;
  ostr << "# columns are min, max, value" << endl;
  for (unsigned int i_hx=0; i_hx<bins; ++i_hx){
    for (unsigned int i_hy=0; i_hy<bins; ++i_hy){
      gsl_histogram2d_get_xrange(h_simple2d, i_hx, &xminval, &xmaxval);
      gsl_histogram2d_get_yrange(h_simple2d, i_hy, &yminval, &ymaxval);
      xbinwidth = xmaxval-xminval;
      ybinwidth = ymaxval-yminval;
      xmeanval=0.5*(xminval+xmaxval);
      ymeanval=0.5*(yminval+ymaxval);
      value2 =sqrt(xmeanval*ymeanval)* xmeanval*ymeanval*gsl_histogram2d_get(h_simple2d, i_hx,i_hy)/(iterations*xbinwidth*ybinwidth);
      ostr2 << xmeanval << " " << ymeanval << " " << value2 << endl;
      if (i_hx==i_hy){
        ostr3<< xminval << " " << xmaxval << " " << value2 << endl;
      }
    }
  }
  ostr2.close();
  ostr3.close();



  system("gnuplot -p 'plotting_simple.p'");
  system("gnuplot -p 'plotting_simple2d.p'");
  system("gnuplot -p 'plotting_simple2d_diag.p'");



/*
  ///Checking bin close to 0.5 and first bin in time
  cout << "Checking fluctuations, takes less than a minute at current number of iterations" <<endl;
  GeneratorInMedium v;
  double cutoff = 1e-10;
  v.set_cutoff(cutoff);
  //unused: double alpha = (4-8*cutoff)/sqrt(cutoff*(1-cutoff));///<integral of simplified kernel
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
        v.generate_branching(1);
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



*/






  return 0;
}

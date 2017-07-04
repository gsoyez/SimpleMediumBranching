#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

#include "event.hh"
#include "generator.hh"
#include "meananderr.hh"

#include <vector>
#include <algorithm>


using namespace std;

int main(){
  cout << "Running"<<endl;

  double tau;
  double epsilon=1e-8;
  double x_min=1e-3;

  GeneratorInMediumSimple g;
  g.set_seed(1);

  GeneratorInMedium g_full;
  g_full.set_seed(1);


  ofstream ostr("energyloss_t3_xmine-3_itere3epsilone-8.dat");
  ostr << "#Time, th average, th variance, simple num average, simple num variance,";
  ostr <<"full num average, full num variance" << endl;

  //for (x_min=0.1; x_min>0.9e-5;x_min/=10){
    //ostr << "#x_min " << x_min<<endl;
    cout << "xmin" << x_min <<endl;

    for (tau=3; tau<3.05;tau+=0.1){
        cout << tau<<endl;
      MeanAndErr energy;
      MeanAndErr energy_full;
      unsigned int iterations = 1000;
      for (unsigned int i = 0; i < iterations; i++ ){
        g.generate_event(tau,epsilon,x_min);

        g_full.generate_event(tau,epsilon,x_min);

        const vector<double> finals=g.event().final_particles();
        double sum = 0;
        vector<double> sortfinals = finals;//not const, can sort
        sort(sortfinals.begin(), sortfinals.end());
        for (unsigned int k = 0; k < sortfinals.size(); k++ ) {
          if (sortfinals[k] > x_min){
            sum += sortfinals[k];
          }
        }
        energy.addEntry(1-sum);


        const vector<double> finals_full=g_full.event().final_particles();
        sum = 0;
        vector<double> sortfinals_full = finals_full;//not const, can sort
        sort(sortfinals_full.begin(), sortfinals_full.end());
        for (unsigned int k = 0; k < sortfinals_full.size(); k++ ) {
          if (sortfinals_full[k] > x_min){
            sum += sortfinals_full[k];
          }
        }
        energy_full.addEntry(1-sum);

        if (i%100 == 0){
          cout << i << endl;
        }
      }


      double erf1=gsl_sf_erf (sqrt(M_PI)*tau);
      double erf2=gsl_sf_erf (2*sqrt(M_PI)*tau);
      double variance = 2*M_PI*tau*(erf1-erf2)+2*exp(-M_PI*tau*tau)
                        - exp(-4*M_PI*tau*tau) - exp(-2*M_PI*tau*tau);

      ostr << tau << "  ";
      ostr << 1-exp(-M_PI*tau*tau) << "  ";
      ostr << variance << "  ";
      ostr << energy.mean()<< "  ";
      ostr << energy.variance() << "  ";
      ostr << energy_full.mean()<< "  ";
      ostr << energy_full.variance() << endl;


    }
  //}



  ostr.close();
  return 0;
}

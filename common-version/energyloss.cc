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



  ofstream ostr("energyloss_survivor.dat");
  ostr << "#Time, th average, th variance, simple num average, simple num variance,";
  ostr <<"full num average, full num variance, simple D2 variance, full D2 variance" << endl;

  ostr << endl << "#x_min " << x_min<<endl;
  cout << "xmin " << x_min <<endl;

  //For tau = 2, seeds 4401 and 6707 gives surviving particle, out of seeds 1-10001
  //One for these, 6707, is also present for tau=3

  x_min=1e-3;
  tau=3;
    cout << "tau " << tau<<endl;
  MeanAndErr energy;

  int seed = 6707;
  g.set_seed(seed);
  g.generate_event(tau,epsilon,x_min);

  const vector<double> finals=g.event().final_particles();
  double sum = 0;
  double xsum = 0;
  double sum2 = 0;
  vector<double> sortfinals = finals;//not const, can sort
  sort(sortfinals.begin(), sortfinals.end());
  for (unsigned int k = 0; k < sortfinals.size(); k++ ) {
    if (sortfinals[k] > x_min){
      sum += sortfinals[k];
      cout << sortfinals[k] << endl;
    }

  }
  if (sum>0){
    cout << "seed " <<seed << ". ";
  }
  energy.addEntry(1-sum);




  ostr.close();
  return 0;
}

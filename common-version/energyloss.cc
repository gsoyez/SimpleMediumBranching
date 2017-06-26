#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

#include <gsl/gsl_math.h>

#include "event.hh"
#include "generator.hh"
#include "meananderr.hh"

#include <vector>
#include <algorithm>


using namespace std;

int main(){
  cout << "Running"<<endl;

  double test_time=0.1;
  double epsilon=1e-8;
  double x_min=1e-4;

  GeneratorInMediumSimple g;
  g.set_seed(1);

  MeanAndErr energy;

  unsigned int iterations = 1000;
  for (unsigned int i = 0; i < iterations; i++ ){
    g.generate_event(test_time,epsilon,x_min);
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
    if (i%10 == 0){
      cout << i << endl;
    }
  }

  double average = energy.mean();
  double error = sqrt(energy.variance());

  cout << "Theory: " << 1-exp(-3.14*test_time*test_time);
  cout << ", average: " << average << ", error*sqrt(3): " << error*sqrt(3) << endl;

  return 0;
}

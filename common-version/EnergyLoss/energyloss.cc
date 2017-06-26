#include <iostream>
#include <fstream>
#include <iomanip>

#include "event.hh"
#include "generator.hh"
#include "meananderr.hh"

using namespace std;

int main(){
  cout << "Running"<<endl;

  double test_time=0.1;
  double epsilon=1e-3;

  GeneratorInMedium gen;
  gen.set_seed(1);
  gen.generate_event(test_time,epsilon);

  return 0;
}

#include <iostream>
#include "event.hh"
#include "generator.hh"
#include <iomanip>
using namespace std;


int main(){

  cout << "Running"<<endl;

  double testTime=0.1;
  double epsilon=1e-3;
	GenerateInMedium gen;
	gen.generateEvent(testTime,epsilon);

	const vector<Particle> testvec=gen.getEvent().getParticles();

  cout <<setw(8) << "index" <<setw(8)<< "parent"<<setw(8) <<"child1"<<setw(8) <<"child2";
  cout <<setw(8) <<"start"<<setw(8) <<"end"<<setw(8) <<"x"<<endl;
  cout << setprecision(2);

  for (unsigned int i = 0; i < testvec.size(); i++ ) {
    cout <<setw(8) << i<<setw(8)<< testvec[i].parent();
    cout <<setw(8) << testvec[i].child1();
    cout <<setw(8) << testvec[i].child2();
    cout <<setw(8) << testvec[i].startTime();
    cout <<setw(8) << testvec[i].endTime();
    cout <<setw(8) << testvec[i].x();
    cout <<endl;
  }

	const vector<Parton> testvec2=gen.getEvent().getPartons();
  cout << endl << "Note that the parton vector below is not filled properly." << endl;
  cout << "We should first decide upon in which order we want to fill it." << endl;
  cout <<setw(8) << "index" <<setw(8)<< "start"<<setw(8) <<"end";
  cout <<setw(8) <<"x"<<setw(8)<< "final"<<endl;
  cout << setprecision(2);

  for (unsigned int i = 0; i < testvec2.size(); i++ ) {
    cout <<setw(8) << i<<setw(8)<< testvec2[i].get_starting_vertex();
    cout <<setw(8) << testvec2[i].get_ending_vertex();
    cout <<setw(8) << testvec2[i].get_x();
    cout <<setw(8) << testvec2[i].is_final();
    cout <<endl;
  }




  return 0;
}

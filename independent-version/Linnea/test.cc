/// std
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime> ///to get some rough comparison between different things
///gsl
#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>

///Own
#include "meananderr.hh"
#include "event.hh"

using namespace std;

int main(){

  cout << "Up and running!"<<endl;

	double testTime=0.2;

	GenerateInMedium gen;
	gen.setSeed(2);
	gen.generateEvent(testTime,1e-3);

	const vector<Particle> testvec=gen.getEvent().getParticles();
  vector<Particle>::const_iterator it2 = testvec.begin();
	cout << "Parent of first with getEvent: " << it2->parent()<<endl;///<causes issues

	cout<<"Total number of particles test : "<< testvec.size()<<endl<<endl;

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
  cout <<setw(8) << "index" <<setw(8)<< "parent"<<setw(8) <<"child1"<<setw(8) <<"child2";
  cout <<setw(8) <<"start"<<setw(8) <<"end"<<setw(8) <<"x"<<endl<<endl;

  int p1index=5;
  int p2index=2;


  ///Test new algorithm for LCA
  int last2;
  clock_t start2 = clock();
  int largest=max(p1index, p2index);
  int smallest=min(p1index, p2index);
  int sibling=largest,parent,current=largest;
  ///Go back in tree, taking advantage of its index structure
  while(sibling>smallest){
    parent=testvec[current].parent();
    if(current==smallest){ ///<One is "along the branch" of the other
      last2 = current;
      break;///<Needed when along branch and sibling larger. Alternative? Break=ugly
    }
    else{
      last2=parent;
    }
    ///Find which child is the sibling
    if (testvec[parent].child1()==current){ ///<child1 same => not the sibling
      sibling=testvec[parent].child2();
    }
    else{ ///They are not the same, so child1 is the sibling
      sibling=testvec[parent].child1();
    }
    current=parent;
  }

  clock_t stop2 = clock();
  double time2 = double(stop2 - start2);
	cout<< "Algorithm 2: last common ancestor of " << p1index << " and " << p2index << " is " << last2;
	cout<<", cpu time " << time2 <<endl;

  return 0;
}

#include "event.hh"
#include <gsl/gsl_math.h>
#include <iostream>
#include <gsl/gsl_integration.h>

using std::cout;
using std::endl;

///Initializing Particle with values that will signify that
///this has no physical meaning. E.g. first particle
///will not have a parent and will get parent = -1.
Particle::Particle(int parent,double startTime,double x,double endTime,
        int child1, int child2){
  setParent(parent);
  setStartTime(startTime);
  setX(x);
  setEndTime(endTime);
  setChild1(child1);
  setChild2(child2);
}
//############################################
///Defining member functions for the class Branching


/// default ctor
Branching::Branching(int init_seed): _t(-1), _z(-1){
  _r = gsl_rng_alloc (gsl_rng_default);
  setSeed(init_seed);
}

/// default dtor
Branching::~Branching(){
  gsl_rng_free(_r);
}

///Using simplified kernel g >= real kernel f in the veto algorithm
/// g = 1/(z(1-z))^(3/2)
/// f = (1-z(1-z))^(5/2)/(z(1-z))^(3/2)
void Branching::generateBranching(double x, double cutoff){
  double z=1,R=1,fgratio=0,t=0;
  cutoff=cutoff/x;
  ///Calculating the integral from cutoff to 1-cutoff
  double alpha = (4-8*cutoff)/sqrt(cutoff*(1-cutoff));

  while (fgratio<R){
    ///Generate t according to g
    t = t - log(gsl_rng_uniform(_r))/(alpha/(2*sqrt(x)));
    ///Generate z according to g. x-dependence of alpha cancels
    double u = gsl_rng_uniform(_r);
    double v = (4*cutoff-2)/sqrt(cutoff*(1-cutoff))+u*alpha*0.5;
    ///z distribution is symmetric around 0.5. Generate below.
    z = 0.5 + v/(2*sqrt(v*v+16));
    ///Decide if rejecting or not
    fgratio = pow((1-z*(1-z)),2.5); ///Rejection factor f/g
    R = gsl_rng_uniform(_r);
  }
  _t=t;
  _z=z;
}


void Branching::setSeed(int new_seed){
	_seed=new_seed;
	gsl_rng_set (_r, _seed);
}

///Defining member functions for the class GenerateInMedium

void GenerateInMedium::generateEvent(double time,double cutoff){
  ///-1 is the default parent, 0 is starting time, x is energy fraction
  _particles.clear();
	Particle firstParticle;
	firstParticle.setStartTime(0);
  firstParticle.setX(1);
  _particles.push_back(firstParticle);
  _endTime=time;
  _cutoff=cutoff;
	_branch();
	_event.setParticles(_particles);
	_event.setCutoff(cutoff);
	_event.setEndTime(time);
}

///Branches the last particle in _particles
void GenerateInMedium::_branch(){
  unsigned int parent=_particles.size()-1;
  double x=_particles[parent].x();///<Relies on pushing back child before branching
  if(x<2*_cutoff){
      return;///<We have reached the bottom of the recursion => we have one particle.
  }
  ///splittingTime, randomly generated
  _v.generateBranching(x, _cutoff);
  double splittingTime = _v.t();
  double z = _v.z();
  double timeLeft=_endTime - _particles[parent].startTime();
  if(timeLeft < splittingTime){
    return;///<We have reached the bottom of the recursion => we have one particle.
  }
  ///No if triggered =>we have time enough for another split and x is large enough
  double currentTime=_particles[parent].startTime() + splittingTime;
  _particles[parent].setEndTime(currentTime);
  ///To add to the vector storing the event
  Particle child1(parent,currentTime,z*x), child2(parent,currentTime,(1-z)*x);
  ///Branching into two
  ///First child:
  _particles.push_back(child1);
  _particles[parent].setChild1(_particles.size()-1);
  _branch();
  ///Second child:
  _particles.push_back(child2);
  _particles[parent].setChild2(_particles.size()-1);
  _branch();
}

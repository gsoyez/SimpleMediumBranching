#include "event.hh"
#include <gsl/gsl_math.h>
#include <iostream>
#include <gsl/gsl_integration.h>
#include "generator.hh"

/// Useful constants

namespace QCD
{
    // Color coefficients
    const double TR = 1./2.;
    const double CA = 3;
    const double CF = 4./3;

    // Coupling constant
    const double alphaS = 0.1;                 ///< Strong coupling constant (taken independent of the energy scale)
    const double alphaSbare = alphaS/(2*M_PI); ///< Normalized strong coupling constant
}

using namespace QCD;



//############################################
// Defining member functions for the class Branching
// which is the class responsible for generating t and z.
// This would be the one for which similar ones should be
// added with different kernels.


/// default ctor
Branching::Branching(int init_seed): _t(-1), _z(-1){
  _r = gsl_rng_alloc (gsl_rng_default);
  setSeed(init_seed);
}

/// default dtor
Branching::~Branching(){
  gsl_rng_free(_r);
}

//------------------------------------------------------------------------
// generate an event
//  - time    maximal time over which we keep branching
//  - cutoff  min x value allowed
//
// Using simplified kernel g >= real kernel f in the veto algorithm
//  g = 1/(z(1-z))^(3/2)
//  f = (1-z(1-z))^(5/2)/(z(1-z))^(3/2)
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
  _partons.clear();
	Particle firstParticle;
  Parton first_parton(-1,0,1);
  _partons.push_back(first_parton);
	firstParticle.setStartTime(0);
  firstParticle.setX(1);
  _particles.push_back(firstParticle);
  _endTime=time;
  _cutoff=cutoff;
	_branch();
	_event.setParticles(_particles);
	_event.setPartons(_partons);
	_event.setVertices(_vertices);
	_event.setCutoff(cutoff);
	_event.setEndTime(time);
}

///Branches the last particle in _particles
void GenerateInMedium::_branch(){
  unsigned int parent=_particles.size()-1;

  unsigned int parent_parton=_partons.size()-1;
  unsigned int parent_vertex=_vertices.size();

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
  ///No "if" triggered =>we have time enough for another split and x is large enough
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

  ///Some structure should be decided on. What is below is just
  ///filler code and does not give anything decent
  Parton first(parent,-1,z*x), second(parent,-1,(1-z)*x);
  _partons.push_back(first);
  _partons.push_back(second);
  Vertex new_vertex(parent_vertex,_partons.size()-2,_partons.size()-1,currentTime);
  _vertices.push_back(new_vertex);

}




///All below here should be integrated into the appropriate classes
const double epsilon = 0.000001;    /// Minimal z value for detectable outgoing gluon
const double time_max = 10;         /// Minimal value for the angle of a resolvable outgoing gluon
const double time_min = 0;          /// Maximal value for the angle, roughly the width of the jet

///                    Defining splitting functions


// Pgg splitting kernel as a function of z.
double Pgg(double z)
{
	return CA*(pow(z,4)+1+pow(1-z,4))/(z*(1-z));
}

// Integration of the Pgg kernel over detectable energy fraction, between epsilon/x=z_min and 1/2.
double IntPgg(double z_min)
{
    // z_min=epsilon/xparent
    return -CA*(11./6.+z_min*(-4+z_min*(1-2./3.*z_min))-4*atanh(1-2*z_min));
}

// Pqqbar splitting kernel as a function of z.
double Pqqbar(double z)
{
    return TR*(pow(z,2)+pow(1-z,2));
}

// Integration of the Pqqbar kernel over detectable energy fraction, between epsilon/x=z_min and 1/2.
double IntPqqbar(double z_min)
{
    return TR*(1./3.+z_min*(-1+(z_min*(1-2./3.*z_min))));
}

// Pqq splitting kernel as a function of z.
double Pqq(double z)
{
    return CF*(1+pow(z,2))/(1-z);
}

// Integration of the Pqqbar kernel between z_min and 1-z_min.
double IntPqq(double z_min)
{
    return CF*(-3./2.+3*z_min+4*atanh(1-2*z_min));
}

///##########################################################################################################################


///                 Random variable generator following the Sudakov probability distribution


// Inverse of the Sudakov probability density
double Inverse_Sudakov(double u, double time_parent, double xparent)
{
    return time_parent-1/(alphaSbare*IntPgg(epsilon/xparent))*log(u);
}

// Random variable generator following the Sudakov probability density
double alea_generator_sudakov(double time_parent, double xparent,gsl_rng *r)
{
	double alea;
	alea = gsl_rng_uniform(r);
	return Inverse_Sudakov(alea,time_parent,xparent);
}


///##########################################################################################################################


///             Random variable generator following the Pgg probability density


// Approximate Pgg density between 1/2 and 1-epsilon for the veto algoritm
double alea_generator_Pgg_approx(double xparent,gsl_rng *r)
{
	double norm, al;
	norm = 2*log(-1+1/(epsilon/xparent)); /// Integral of 2/(z(1-z)) between epsilon/x and 1/2
	al = gsl_rng_uniform(r)*norm;
	return 1/(1+exp(-al/2.)); /// Generates a number between 1/2 and 1-epsilon/x
}

// Implementation of the veto algorithm, Pgg density between 1/2 and 1-epsilon
double alea_generator_Pgg(double xparent,gsl_rng *r)
{
	double output, z_candidate, alea;
	while(1)
    {
		z_candidate = alea_generator_Pgg_approx(xparent,r);
		alea = gsl_rng_uniform(r)*CA*2/(z_candidate*(1-z_candidate));
		output = z_candidate;
		if(alea<=Pgg(z_candidate))
        {
			break;
        }
    }
	return output; /// Generates a number between 1/2 and 1-epsilon/
}



///##########################################################################################################################


///              Random variable for the medium induced K(z) probability density


// Approximate K density between epsilon and 1/2
double alea_generator_K_approx(double xparent, gsl_rng *r)
{
    double z_min = epsilon/xparent;
    double norm;
    norm = (2*(1-2*z_min))/sqrt(z_min*(1-z_min)); /// Integral of 1/(z*(1-z)0^(3/2) between epsilon/x and 1/2
    double al = gsl_rng_uniform(r)*norm;
    al*=al;
    return (16+al-sqrt(16*al+al*al))/(2*(16+al)); /// Generates a number between epsilon/x and 1/2
}


///##########################################################################################################################







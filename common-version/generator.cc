#include "event.hh"
#include <gsl/gsl_math.h>
#include <iostream>
#include <cassert>
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
using namespace std;

//Defining member functions for the class GeneratorInMedium

//------------------------------------------------------------------------
// implementation of the base class methods
//------------------------------------------------------------------------
// default ctor
GeneratorBase::GeneratorBase(int init_seed): _t(-1), _z(-1){
  _r = gsl_rng_alloc (gsl_rng_default);
  set_seed(init_seed);

  _recursive_branching = true;
}

// default dtor
GeneratorBase::~GeneratorBase(){
  gsl_rng_free(_r);
}

// set the random seed
void GeneratorBase::set_seed(int new_seed){
  _seed=new_seed;
  gsl_rng_set (_r, _seed);
}

// generate an event
//  - time    maximal time over which we keep branching
//  - cutoff  min x value allowed for emissions
//  - xmin    min x value allowed for further branching
//            (-ve means 2*cutoff i.e. no effect)
void GeneratorBase::generate_event(double time,double cutoff,double xmin){
  //--------------------------------------------------
  // initialise the event and local variables
  
  // -1 is the default parent, 0 is starting time, x is energy fraction
  _particles.clear();
  Particle first_particle;
  first_particle.set_start_time(0);
  first_particle.set_x(1);
  _particles.push_back(first_particle);
  _end_time=time;
  _cutoff=cutoff;
  set_xmin(xmin);

  // start the branching
  if (_recursive_branching){
    _branch_last_recursively();
  } else {
    _branch_non_recursively();
  }    

  // create the final event
  _event.set_particles(_particles);
  _event.set_cutoff(cutoff);
  _event.set_end_time(time);
}

// Branches the last particle in _particles, then iterate
void GeneratorBase::_branch_last_recursively(){
  //GS comments: this can be largely simplified:
  // . splitting_time and _z are not needed (use _t and _z)
  // . time_left is not needed and should be replaced by current_time

  unsigned int parent=_particles.size()-1;

  double x=_particles[parent].x();///< Relies on pushing back child before branching
  if(x<_xmin){
    return;///< We have reached the bottom of the recursion => we have one particle.
  }
  /// SplittingTime, randomly generated
  generate_branching(x);
  double splitting_time=_t;
  double z=_z;

  double time_left=_end_time - _particles[parent].start_time();
  if(time_left < splitting_time){
    return;///< We have reached the bottom of the recursion => we have one particle.
  }
  /// No "if" triggered =>we have time enough for another split and x is large enough
  double current_time=_particles[parent].start_time() + splitting_time;
  _particles[parent].set_end_time(current_time);
  /// To add to the vector storing the event
  Particle child1(parent,current_time,z*x), child2(parent,current_time,(1-z)*x);

  /// Branching into two
  /// First child:
  _particles.push_back(child1);
  _particles[parent].set_child1(_particles.size()-1);
  _branch_last_recursively();
  ///Second child:
  _particles.push_back(child2);
  _particles[parent].set_child2(_particles.size()-1);
  _branch_last_recursively();
}


// generate a whole cascade stating form the existing particles in _particle
void GeneratorBase::_branch_non_recursively(){
  // we start with a vector that will hold the indices of the
  // particles that ned to be forther branched.
  // Initially, it should contain all the particles in the event
  vector<unsigned int> indices_to_branch(_particles.size());
  for (unsigned int i=0; i<_particles.size(); ++i){
    if (_particles[i].x() > _xmin) indices_to_branch[i] = i;
  }

  // not branch until we've exhausted the list of particles to branch
  unsigned int current_position = 0;
  while (current_position < indices_to_branch.size()){
    // get the particle that we're supposed to branch
    unsigned int particle_index = indices_to_branch[current_position];
    const Particle & to_branch = _particles[particle_index];
    double x = to_branch.x();

    // generate (randomly) teh splitting time and momentum fraction
    generate_branching(x);

    // it we've exceeded time, discard the branching
    double branching_time = to_branch.start_time() + _t;
    if (branching_time > _end_time){
      // just go to ythe next branching particle
      ++current_position;
      continue;
    }

    // generate the new particles
    double x1 = _z*x;
    double x2 = (1-_z)*x;
    Particle child1(particle_index, branching_time, x1);
    Particle child2(particle_index, branching_time, x2);

    // update the event record
    _particles.push_back(child1);
    _particles.push_back(child2);
    _particles[particle_index].set_end_time(branching_time);
    _particles[particle_index].set_child1(_particles.size()-2);
    _particles[particle_index].set_child2(_particles.size()-1);

    // proceed w next branchings
    if (x1>_xmin){
      // we'll deal w the branching of x1 right now
      indices_to_branch[current_position] = _particles.size()-2;

      if (x2>_xmin) // append it for future branching
        indices_to_branch.push_back(_particles.size()-1);
    } else {
      if (x2>_xmin) // deal with it immediately
        indices_to_branch[current_position] = _particles.size()-1;
      else // non of the 2 new particles need further branching
        ++current_position;
    }     
  }    
}

//------------------------------------------------------------------------
// implementation of the in-medium branching for the full kernel
//------------------------------------------------------------------------

//------------------------------------------------------------------------
// generate an event 
//  - time    maximal time over which we keep branching
//  - cutoff  min x value allowed
//
// Using simplified kernel g >= real kernel f in the veto algorithm
//  g = 1/(z(1-z))^(3/2)
//  f = (1-z(1-z))^(5/2)/(z(1-z))^(3/2)
//------------------------------------------------------------------------
void GeneratorInMedium::generate_branching(double x){
  double z=1,R=1,fgratio=0,t=0;
  double cutoff=_cutoff/x;
  double root_x=sqrt(x);
  /// Calculating the integral from cutoff to 1-cutoff
  double alpha = (2-4*cutoff)/sqrt(cutoff*(1-cutoff));

  while (fgratio<R){
    /// Generate t according to g
    t = t - log(gsl_rng_uniform(_r))/(alpha/(root_x));
    /// Generate z according to g. x-dependence of alpha cancels
    double u = gsl_rng_uniform(_r);
    double v = (u-1)*alpha;
    /// z distribution is symmetric around 0.5. Generate below.
    z = 0.5 + v/(2*sqrt(v*v+16));
    /// Decide if rejecting or not
    fgratio = pow((1-z*(1-z)),5); /// Rejection factor f/g, squared
    u = gsl_rng_uniform(_r);
    R=u*u;
  }
  _t=t;
  _z=z;
}


//------------------------------------------------------------------------
// implementation of the in-medium branching for the simplified kernel
//------------------------------------------------------------------------

void GeneratorInMediumSimple::generate_branching(double x){
  double z, t; //=1,R=1,fgratio=0,t=0;
  double cutoff=_cutoff/x;
  double root_x=sqrt(x);
  /// Calculating the integral from cutoff to 1-cutoff
  double alpha = (2-4*cutoff)/sqrt(cutoff*(1-cutoff));

  // Generate t according to g
  //t = t - log(gsl_rng_uniform(_r))/(alpha/(root_x));
  t = - log(gsl_rng_uniform(_r))/(alpha/(root_x));
  // Generate z according to g. x-dependence of alpha cancels
  double u = gsl_rng_uniform(_r);
  double v = (u-1)*alpha;
  // z distribution is symmetric around 0.5. Generate below.
  z = 0.5 + v/(2*sqrt(v*v+16));

  _t=t;
  _z=z;
}


//-----------------------------------------------------------------------------------------//




/// All below here should be integrated into the appropriate classes
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

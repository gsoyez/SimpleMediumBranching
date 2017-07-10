#ifndef GENERATOR_HH_INCLUDED
#define GENERATOR_HH_INCLUDED

#include "event.hh"

//---------------------------------------------------------------------------------------------//
//
// In this file we should have all classes that have to do with generators
//
//---------------------------------------------------------------------------------------------//


//Paul's Pgg Pgq etc functions should go here into a class. Use a class for generators and
//then derived classes

/// \class GeneratorBase
/// Base class for generating events. We'll derive GeneratorInMedium
/// and GeneratorInMedium from this.
class GeneratorBase{
public:
  /// default ctor
  GeneratorBase(int init_seed=1);

  /// default dtor
  virtual ~GeneratorBase();

  /// generate an event
  ///  \param time    maximal time over which we keep branching
  ///  \param cutoff  min x value allowed for emissions
  ///  \param xmin    min x value allowed for further branching
  ///                 (-ve means 2*cutoff i.e. no effect)
  void generate_event(double time, double cutoff,double xmin=-1.0);

  /// Return the generated event
  const Event& event() const {return _event;}

  /// related to branching
  virtual void generate_branching(double x)=0;
  double t() const {return _t;}
  double z() const {return _z;}

  /// get/set the random seed
  int seed() const { return _seed;}
  void set_seed(int new_seed);

  /// Set and get the epsilon used for the event
  void set_cutoff(double value) {_cutoff = value;}
  double cutoff() const {return _cutoff;}

  /// allow to stop branching at a given x value (by default set to 2
  /// epsilon). If the value is -ve, reset to the default
  void set_xmin(double value){
    _xmin = (value<=0) ? 2*_cutoff : value;
  }
  double xmin() const { return _xmin; }

  /// switchies between recursive and non-recursive implementation of
  /// the branching
  ///
  /// The recursive strategy is a bit more elegant in terms of coding
  /// and has the advantage of generating trees organised so that it
  /// is easy to extract subtrees
  ///
  /// Conversely, the non-recursve strategy avoid potential stack
  /// overflow problems when going to large times os small x.
  void set_recursive_branching(bool value){ _recursive_branching = value;}
  bool recursive_branching() const { return _recursive_branching;}
  
protected:
  /// stores the last generated event
  Event _event;

  /// stores the list of particles in the event which is currently being generated
  std::vector<Particle> _particles;

  /// branch the last particle in the "_particles" vector, then
  /// recursively branch the result.
  void _branch_last_recursively();

  /// starting from the particles in the _particle vector, branch
  /// (non-recursively) until we've reached teh maximal time.
  void _branch_non_recursively();
  
  // the parameters of the generator
  double _cutoff;   ///< min x value at which we can emit particles
  double _xmin;;    ///< min x value for further branching
  double _end_time; ///< max time up to which one should branch 

  // Related to branching
  gsl_rng *_r; ///< random number generator (GSL)
  int _seed;   ///< random seed
  double _t;   ///< time at which the last branging occured
  double _z;   ///< z fraction for the last branching

  // options to branch recursively or not
  bool _recursive_branching;
};


/// \class GeneratorInMedium
/// Generate an event according to the in-medium evolution from
/// .... This version implementes the full evolution kernel.
class GeneratorInMedium : public GeneratorBase{
public:
  /// default ctor
  GeneratorInMedium(int init_seed=1) : GeneratorBase(init_seed){}

  /// default dtor
  virtual ~GeneratorInMedium(){}

  /// we need to overload the generation of a single branching
  virtual void generate_branching(double x);
};



/// \class GeneratorInMediumSimple
/// Generate an event according to the in-medium evolution from
/// .... This version implementes tyhe simplified evolution kernel.
///
class GeneratorInMediumSimple : public GeneratorBase{
public:
  /// default ctor
  GeneratorInMediumSimple(int init_seed=1) : GeneratorBase(init_seed){}

  /// default dtor
  virtual ~GeneratorInMediumSimple(){}

  /// we need to overload the generation of a single branching
  virtual void generate_branching(double x);
};


#endif // GENERATOR_HH_INCLUDED

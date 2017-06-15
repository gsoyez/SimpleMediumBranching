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


/// \class GeneratorInMedium
/// Generate an event. Returns a member of the class Event.
class GeneratorInMedium{
public:
  /// default ctor
  GeneratorInMedium(int init_seed=1);

  /// default dtor
  ~GeneratorInMedium();

  /// generate an event
  ///  \param time    maximal time over which we keep branching
  ///  \param cutoff  min x value allowed for emissions
  ///  \param xmin    min x value allowed for further branching
  ///                 (-ve means 2*cutoff i.e. no effect)
  void generate_event(double time, double cutoff,double xmin=-1.0);

  /// Return the generated event
  const Event& event() const {return _event;}

  /// related to branching
  void generate_branching(double x);
  double t() const {return _t;}
  double z() const {return _z;}

  /// get the random seed
  int seed() const { return _seed;}
  /// set the random seed
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
  
private:
  Event _event;
  std::vector<Particle> _particles;
  void _branch();
  double _cutoff; ///< sets the min x value at which we can emit particles
  double _xmin;;  ///< sets the min x value for further branching
  double _end_time;

  /// Related to branching
  gsl_rng *_r; ///< random number generator (GSL)
  int _seed;   ///< random seed
  double _t;
  double _z;
};



/// \class GeneratorInMediumSimple
/// Generate an event. Returns a member of the class Event.
class GeneratorInMediumSimple{
public:
  /// default ctor
  GeneratorInMediumSimple(int init_seed=1);

  /// default dtor
  ~GeneratorInMediumSimple();

  /// generate an event
  ///  \param time    maximal time over which we keep branching
  ///  \param cutoff  min x value allowed for emissions
  ///  \param xmin    min x value allowed for further branching
  ///                 (-ve means 2*cutoff i.e. no effect)
  void generate_event(double time, double cutoff, double xmin=-1.0);

  /// Return the generated event
  const Event& event() const {return _event;}

  /// related to branching
  void generate_branching(double x);
  double t() const {return _t;}
  double z() const {return _z;}
  /// get the random seed
  int seed() const { return _seed;}
  /// set the random seed
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

private:
  Event _event;
  std::vector<Particle> _particles;
  void _branch();
  double _cutoff; ///< sets the min x value at which we can emit particles
  double _xmin;;  ///< sets the min x value for further branching
  double _end_time;

  /// Related to branching
  gsl_rng *_r; ///< random number generator (GSL)
  int _seed;   ///< random seed
  double _t;
  double _z;
};


#endif // GENERATOR_HH_INCLUDED

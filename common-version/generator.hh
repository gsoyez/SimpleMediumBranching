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
  ///  \param cutoff  min x value allowed
  void generate_event(double time, double cutoff);

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

private:
  Event _event;
  std::vector<Particle> _particles;
  void _branch();
  double _cutoff;
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
  ///  \param cutoff  min x value allowed
  void generate_event(double time, double cutoff);

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

private:
  Event _event;
  std::vector<Particle> _particles;
  void _branch();
  double _cutoff;
  double _end_time;

  /// Related to branching
  gsl_rng *_r; ///< random number generator (GSL)
  int _seed;   ///< random seed
  double _t;
  double _z;
};


#endif // GENERATOR_HH_INCLUDED

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


/// \class GenerateInMedium
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
  void generateEvent(double time, double cutoff);

  const Event& event() const {return _event;}///< Return the generated event

  /// related to branching
  void generateBranching(double x, double cutoff);
  double t() const {return _t;}
  double z() const {return _z;}
  /// get the random seed
  int seed() const { return _seed;}
  /// set the random seed
  void setSeed(int new_seed);


private:
  Event _event;
  std::vector<Particle> _particles;
  void _branch();
  double _cutoff;
  double _endTime;

  /// Related to branching
  gsl_rng *_r; ///< random number generator (GSL)
  int _seed;   ///< random seed
  double _t;
  double _z;
};


#endif // GENERATOR_HH_INCLUDED

#ifndef GENERATOR_HH_INCLUDED
#define GENERATOR_HH_INCLUDED

#include "event.hh"
//---------------------------------------------------------------------------------------------//

///In this file we should have all classes that have to do with generators

//---------------------------------------------------------------------------------------------//


//Paul's Pgg Pgq etc functions should go here into a class. Use a class for generators and
//then derived classes



/// \class Branching
/// Generate random branching times and z values
/// The class uses GSL for random numbers.
class Branching{
public:
  /// default ctor
  Branching(int init_seed=1);
  /// default dtor
  ~Branching();
  void generateBranching(double x, double cutoff);
  double t() const {return _t;}
  double z() const {return _z;}
  /// get the random seed
  int seed() const { return _seed;}
  /// set the random seed
  void setSeed(int new_seed);

private:
  gsl_rng *_r; ///< random number generator (GSL)
  int _seed;   ///< random seed
  double _t;
  double _z;
};



///\class GenerateInMedium
///Generate an event. Returns a member of the class Event.
class GenerateInMedium{
public:
  void generateEvent(double time, double cutoff);
  const Event getEvent() const {return _event;}///<Return the generated event:
	void setSeed(int newSeed){_v.setSeed(newSeed);}
  int getSeed(){return _v.seed();}

private:
  Event _event;
  std::vector<Particle> _particles;
  std::vector<Parton> _partons;
  std::vector<Vertex> _vertices;
  void _branch();
  Branching _v;
  double _cutoff;
  double _endTime;
};


#endif // GENERATOR_HH_INCLUDED

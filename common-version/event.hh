#ifndef EVENT_HH_INCLUDED
#define EVENT_HH_INCLUDED

#include <vector>
#include "gsl/gsl_rng.h"


//---------------------------------------------------------------------------------------------//
/// In this file we should have all classes that have to do with events
//---------------------------------------------------------------------------------------------//



class Particle{
public:
  Particle(int parent=-1,double startTime=-1,double x=-1,double endTime=-1,
           int child1=-1, int child2=-1);

  /// Set and get parent particle
  void setParent(int value) { _parent=value;}
  int parent() const { return _parent;}

  /// Set and get time when particle is created through branching
  void setStartTime(double value) { _startTime=value;}
  double startTime() const { return _startTime;}

  /// Set and get x(energy fraction compared to leading particle)
  void setX(double value) { _x=value;}
  double x() const { return _x;}

  /// Set and get time when particle branches
  void setEndTime(double value) { _endTime=value;}
  double endTime() const { return _endTime;}

  /// Set and get the child particles the particle branches into
  void setChild1(double value) { _child1=value;}
  void setChild2(double value) { _child2=value;}
  int child1() const { return _child1;}
  int child2() const { return _child2;}

  /// returns true if a particle is final
  bool is_final() const{ return _child1==-1;}

private:
  /// All these are initialized to -1
  int _parent; ///< Not unsigned, want -1 to be an option for the first one
  double _startTime;
  double _x; ///<Fraction of energy compared to leading particle
  double _endTime;
  int _child1;
  int _child2;
};


/// \class Event
/// Generate an event that has branched for a
/// given time, which is given as input
class Event{
public:
  /// Set and get the epsilon used for the event
  void set_cutoff(double value) {_cutoff = value;}
  double cutoff() const {return _cutoff;}

  /// Set and get the time during which the leading particle has branched
  void set_end_time(double value) {_endTime = value;}
  double end_time() const {return _endTime;}
  unsigned int get_final_number() const {return (_particles.size()+1)/2;}///<Number of final particles

  /// Set and get the vector storing all (intermediate and final) particles
  /// When returning vector the compiler should optimize so it shouldn't be an issue
  void set_particles(std::vector<Particle> value){_particles = value;}
  const std::vector<Particle>& particles() const{return _particles;}

  /// Function that returns list of x ie energy fractions for the final partons
  const std::vector<double> final_particles() const;

private:
  std::vector<Particle> _particles;
  double _cutoff;
  double _endTime;
};

#endif // EVENT_HH_INCLUDED

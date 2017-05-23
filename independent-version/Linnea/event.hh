#ifndef EVENT_HH_INCLUDED
#define EVENT_HH_INCLUDED

#include <vector>
#include "gsl/gsl_rng.h"

class Particle{
public:
    Particle(int parent=-1,double startTime=-1,double x=-1,double endTime=-1,
        int child1=-1, int child2=-1);
    ///Set and get parent particle
    void setParent(int value) { _parent=value;}
    int parent() const { return _parent;}
    ///Set and get time when particle is created through branching
    void setStartTime(double value) { _startTime=value;}
    double startTime() const { return _startTime;}
    ///Set and get x(energy fraction compared to leading particle)
    void setX(double value) { _x=value;}
    double x() const { return _x;}
    ///Set and get time when particle branches
    void setEndTime(double value) { _endTime=value;}
    double endTime() const { return _endTime;}
    ///Set and get the child particles the particle branches into
    void setChild1(double value) { _child1=value;}
    void setChild2(double value) { _child2=value;}
    int child1() const { return _child1;}
    int child2() const { return _child2;}

private:
    ///All these are initialized to -1
    int _parent; ///< Not unsigned, want -1 to be an option for the first one
    double _startTime;
    double _x; ///<Fraction of energy compared to leading particle
    double _endTime;
    int _child1;
    int _child2;
};


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


/// \class Event
///Generate an event that has branched for a
///given time, which is given as input
class Event{
public:
    Event();
    void setCutoff(double cutoff);
    double getCutoff() const {return _cutoff;}
    void generateEvent(double time);///nope
	unsigned int getN() const;
	const std::vector<Particle>& getParticles() const;
	///need setParticles
	void setSeed(int newSeed){_v.setSeed(newSeed);}///nope
    int getSeed(){return _v.seed();}///nope
    std::vector<Particle> particles;///<public for now, consider options
    const std::vector<Particle>& testParticles() const {return particles;}///temporary test
private:
    std::vector<Particle> _particles;
    void _branch();///nope
    Branching _v;///nope
    double _cutoff;
    double _endTime;
};



///\class GenerateInMedium
///Generate an event. Returns a member of the class Event.
class GenerateInMedium{
public:
    GenerateInMedium();
    ///Return the generated event:
    void generateEvent(double time, double cutoff);
    const Event event() const {return _event;}

	void setSeed(int newSeed){_v.setSeed(newSeed);}
    int getSeed(){return _v.seed();}
    Event _event; ///should not be public but want to check things
private:

    void _branch();
    Branching _v;
    double _cutoff;
    double _endTime;
};







#endif // EVENT_HH_INCLUDED

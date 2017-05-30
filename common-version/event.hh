#ifndef EVENT_HH_INCLUDED
#define EVENT_HH_INCLUDED

#include <vector>
#include "gsl/gsl_rng.h"


//---------------------------------------------------------------------------------------------//

///In this file we should have all classes that have to do with events

//---------------------------------------------------------------------------------------------//

class Vertex
{
  public:
  // Constructors

  /// default ctor
  /// put everything to -1 to mark them as uninitialised
  Vertex()
  {this->index_parent=-1; this->index_first=-1; this->index_second=-1; this->time=0;};

  /// ctor w initialisation
  Vertex(int index_parent, int index_first, int index_second, double time)
  {this->index_parent=index_parent; this->index_first=index_first; this->index_second=index_second; this->time=time;}

  // Getters
  int get_index_parent() const { return this->index_parent;}
  int get_index_first() const { return this->index_first;}
  int get_index_second() const { return this->index_second;}
  double get_time() const { return this->time;}

  // Setters
  void set_index_parent(int index_parent) {this->index_parent=index_parent;}
  void set_index_first(int index_first) {this->index_first=index_first;}
  void set_index_second(int index_second) {this->index_second=index_second;}
  void set_time(double time) {this->time=time;}

	private:

  int index_parent;  ///< Index in a given list of the parent parton.
	int index_first;   ///< Index in a given list of the first parton going out the vertex
	int index_second;  ///< Index in a given list of the second parton going out the vertex
	double time;       ///< Square of the angle associated with the vertex
};

class Parton
{
  public:

  // Constructors
  Parton()
  {this->starting_vertex=-1;this->ending_vertex=-1; this->x=1;};
  Parton(int starting_vertex,int ending_vertex, double x)
  {this->starting_vertex=starting_vertex; this->ending_vertex=ending_vertex; this->x=x;};

  // Getters
  int get_starting_vertex() const { return this->starting_vertex;}
  double get_ending_vertex() const { return this->ending_vertex;}
  double get_x() const {return this-> x;}

  // Setters
  void set_starting_vertex(int starting_vertex) {this->starting_vertex=starting_vertex;}
  void set_ending_vertex(int starting_vertex) {this->starting_vertex=starting_vertex;}
  void set_x(double x) {this->x=x;}

  /// returns true if a particle is final
  bool is_final() const{ return ending_vertex==-1; }

	private:

	int starting_vertex;   ///< Index in given list of the vertex which created the parton
	double ending_vertex;  ///< Index in given list of the vertex which deleted the parton
	double x;              ///< Total fraction of energy x carried by the parton
};


//---------------------------------------------------------------------------------------------//

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
  /// returns true if a particle is final
  bool is_final() const{ return _child1==-1;}

private:
  ///All these are initialized to -1
  int _parent; ///< Not unsigned, want -1 to be an option for the first one
  double _startTime;
  double _x; ///<Fraction of energy compared to leading particle
  double _endTime;
  int _child1;
  int _child2;
};


/// \class Event
///Generate an event that has branched for a
///given time, which is given as input
class Event{
public:
  ///Set and get the epsilon used for the event
  void setCutoff(double value) {_cutoff = value;}
  double getCutoff() const {return _cutoff;}
  ///Set and get the time during which the leading particle has branched
  void setEndTime(double value) {_endTime = value;}
  double getEndTime() const {return _endTime;}
	unsigned int getN() const {return (_particles.size()+1)/2;}///<Number of final particles
	///Set and get the vector storing all (intermediate and final) particles
	///When returning vectors the compiler should optimize so it shouldn't be an issue
	void setParticles(std::vector<Particle> value){_particles = value;}
	const std::vector<Particle> getParticles() const{return _particles;}
	void setPartons(std::vector<Parton> value){_partons = value;}
	const std::vector<Parton> getPartons() const{return _partons;}
	void setVertices(std::vector<Vertex> value){_vertices = value;}
	const std::vector<Vertex> getVertices() const{return _vertices;}
	///Function that returns list of x ie energy fractions for the final partons
  const std::vector<double> get_final() const;
private:
  std::vector<Particle> _particles;
  std::vector<Parton> _partons;
  std::vector<Vertex> _vertices;
  double _cutoff;
  double _endTime;
};

#endif // EVENT_HH_INCLUDED

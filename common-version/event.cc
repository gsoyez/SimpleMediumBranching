#include "event.hh"

///Initializing Particle with values that will signify that
///this has no physical meaning. E.g. first particle
///will not have a parent and will get parent = -1.
Particle::Particle(int parent,double startTime,double x,double endTime,
        int child1, int child2){
  setParent(parent);
  setStartTime(startTime);
  setX(x);
  setEndTime(endTime);
  setChild1(child1);
  setChild2(child2);
}


///Member function of Event that returns final partons' x-values
const std::vector<double> Event::final_particles() const{
  std::vector<double> x_values;
  for (unsigned int i = 0; i < _particles.size(); i++ ){
    if (_particles[i].is_final()){
      x_values.push_back(_particles[i].x());
    }
  }
  return x_values;
}

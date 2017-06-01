#include "event.hh"

///Initializing Particle with values that will signify that
///this has no physical meaning. E.g. first particle
///will not have a parent and will get parent = -1.
Particle::Particle(int parent,double start_time,double x,double end_time,
        int child1, int child2){
  set_parent(parent);
  set_start_time(start_time);
  set_x(x);
  set_end_time(end_time);
  set_child1(child1);
  set_child2(child2);
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

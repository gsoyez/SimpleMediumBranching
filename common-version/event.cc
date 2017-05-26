#include "event.hh"
#include <gsl/gsl_math.h>
#include <iostream>
#include <gsl/gsl_integration.h>


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

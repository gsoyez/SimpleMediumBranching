#include "meananderr.hh"
#include <gsl/gsl_math.h>

///initializing the values
MeanAndErr::MeanAndErr(){
  _counter = 0;
  _value=0;
  _valueSquared=0;
}


///Updating values
void MeanAndErr::addEntry(double myValue){
  _counter +=1;
  _value += myValue;
  _valueSquared += myValue*myValue;
}


///Functions to get the relevant values
double MeanAndErr::mean() const{
    return (double)_value/_counter;
}

double MeanAndErr::variance() const{
    double n = (double)_value/_counter;
    double n2 = (double)_valueSquared/_counter;
    return (double)_counter/(_counter-1)*(n2 - n*n);
}

double MeanAndErr::error() const{
    double n = (double)_value/_counter;
    double n2 = (double)_valueSquared/_counter;
    return sqrt((double)1/(_counter-1)*(n2 - n*n));
}


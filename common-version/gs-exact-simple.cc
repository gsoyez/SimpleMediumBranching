#include <gsl/gsl_integration.h>
#include <cmath>
#include <CmdLine.hh>
#include <iostream>
#include <fstream>

double tmax;

//========================================================================
double D(double x, void *unused){
  return tmax/sqrt(x*(1-x)*(1-x)*(1-x))*exp((-M_PI*tmax*tmax)/(1-x));
}

int main(int argc, char*argv[]){
  CmdLine cmd(argc, argv);
  tmax = cmd.value("-tmax", 1.0);
  double xmin = cmd.value("-xmin", 0.0001);
  unsigned int nbin = cmd.value<unsigned int>("-nbin", 80);
  string out = cmd.value<string>("-out");

  //----------------------------------------------------------------------
  // a header with or setup
  ofstream ostr(out.c_str());
  ostr << "# Ran: " << cmd.command_line() << endl;
  ostr << "#" << endl;
  ostr << "# tmax = " << tmax << endl;
  ostr << "# xmin = " << xmin << endl;

  //----------------------------------------------------------------------
  // setup calculation tools
  unsigned int wsize = 100000;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(wsize);
  gsl_function f;
  f.function = &D;

  //----------------------------------------------------------------------
  // loop over the bins in logx
  double lxmin = log(1.0/xmin);
  for (unsigned int i=0;i<nbin; ++i){
    double l1 = exp(-(i+0.0)*lxmin/nbin);
    double l2 = exp(-(i+1.0)*lxmin/nbin);

    double res, err;
    gsl_integration_qags(&f, l1, l2, 1.0e-8, 1.0e-8, wsize, w, &res, &err);
    ostr << l1 << " " << 0.5*(l1+l2) << " " << l2 << " " 
         << res/(l2-l1) << " " << err/(l2-l1) << endl;
  }
  

  gsl_integration_workspace_free(w);

  return 0;
}


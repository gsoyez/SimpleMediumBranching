#include <gsl/gsl_integration.h>
#include <cmath>
#include <CmdLine.hh>
#include <iostream>
#include <fstream>

double t;

//========================================================================
double D(double x, void *unused){
  return t/sqrt(x*(1-x)*(1-x)*(1-x))*exp((-M_PI*t*t)/(1-x));
}

int main(int argc, char*argv[]){
  CmdLine cmd(argc, argv);
  double tmax = cmd.value("-tmax", 1.0);
  unsigned int nt = cmd.value<unsigned int>("-nt", 10);
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

  for (unsigned int it=0; it<nt; ++it){
    t = (it+1.0)/nt*tmax;
    ostr << "# D(x,t=" << t << ")" << endl;

    for (unsigned int i=0;i<nbin; ++i){
      double l1 = (i+0.0)*lxmin/nbin;
      double l2 = (i+1.0)*lxmin/nbin;

      double x1 = exp(-l1);
      double x2 = exp(-l2);

      double res, err;
      gsl_integration_qags(&f, x2, x1, 1.0e-8, 1.0e-8, wsize, w, &res, &err);
      ostr << l1 << " " << 0.5*(l1+l2) << " " << l2 << " " 
           << res/(x1-x2) << " " << err/(x1-x2) << endl;
    }

    ostr << endl << endl;
  }
  

  gsl_integration_workspace_free(w);

  return 0;
}


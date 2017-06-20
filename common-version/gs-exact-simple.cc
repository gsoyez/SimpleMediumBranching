#include <gsl/gsl_integration.h>
#include <cmath>
#include <CmdLine.hh>
#include <iostream>
#include <fstream>

double t;

// helpers for GSL (avoid reallocating the workspace every time
const unsigned int wsize = 100000;
gsl_integration_workspace *wnest;
gsl_function fnest;

//========================================================================
double D(double x, void *unused){
  return t/sqrt(x*(1-x)*(1-x)*(1-x))*exp((-M_PI*t*t)/(1-x));
}

double D2(double x, void *xpptr){
  double xp = (double) (* (double*)xpptr);
  if (x+xp>=1.0) return 0.0;
  double tmp = (-M_PI*t*t)/(1-x-xp);
  return 1.0/(2*M_PI*sqrt(x*xp*(1-x-xp)))*(exp(tmp)-exp(4*tmp));
}

struct XBounds{
  double xmin, xmax;
};
  
double D2_first(double xp, void *xbounds_ptr){
  const XBounds & xbounds = (XBounds) (* (XBounds*)xbounds_ptr);
  fnest.params = (void*) (&xp);
  double res, err;
  gsl_integration_qags(&fnest, xbounds.xmin, xbounds.xmax, 1.0e-8, 1.0e-8, wsize, wnest, &res, &err);
  return res;
}


int main(int argc, char*argv[]){
  CmdLine cmd(argc, argv);
  double tmax = cmd.value("-tmax", 1.0);
  unsigned int nt = cmd.value<unsigned int>("-nt", 10);
  double xmin = cmd.value("-xmin", 0.0001);
  unsigned int nbin = cmd.value<unsigned int>("-nbin", 80);
  unsigned int nD2  = cmd.value<unsigned int>("-nD2", 40);
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
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(wsize);
  gsl_function f;
  f.function = &D;

  wnest = gsl_integration_workspace_alloc(wsize);
  fnest.function = &D2;

  //----------------------------------------------------------------------
  // D(x,t)
  // loop over the t values and the bins in logx
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

  //----------------------------------------------------------------------
  // D2(x,x;t=1)
  // fix t=tmax and loop over the bins in logx
  f.function = &D2_first;
  XBounds xbounds;
  f.params = (void*) (&xbounds);

  t = tmax;
  ostr << "# D2(x,t=" << t << ")" << endl;

  for (unsigned int i=0;i<nbin; ++i){
    double l1 = (i+0.0)*lxmin/nbin;
    double l2 = (i+1.0)*lxmin/nbin;
    
    double x1 = exp(-l1);
    double x2 = exp(-l2);
    
    xbounds.xmin = x2;
    xbounds.xmax = x1;
    
    double res, err;
    gsl_integration_qags(&f, x2, x1, 1.0e-8, 1.0e-8, wsize, w, &res, &err);
    ostr << l1 << " " << 0.5*(l1+l2) << " " << l2 << " " 
         << res/((x1-x2)*(x1-x2)) << " " << err/((x1-x2)*(x1-x2)) << endl;
  }
  
  ostr << endl << endl;

  //----------------------------------------------------------------------
  // D2(x,x';t=1)
  // fix t=tmax and loop over the bins in logx
  ostr << "# D2full(x,x',t=" << t << ")" << endl;

  for (unsigned int i=0;i<nD2; ++i){
    double l1 = (i+0.0)*lxmin/nD2;
    double l2 = (i+1.0)*lxmin/nD2;
    
    double x1 = exp(-l1);
    double x2 = exp(-l2);
    
    xbounds.xmin = x2;
    xbounds.xmax = x1;

    for (unsigned int ip=0;ip<nD2; ++ip){
      double lp1 = (ip+0.0)*lxmin/nD2;
      double lp2 = (ip+1.0)*lxmin/nD2;
    
      double xp1 = exp(-lp1);
      double xp2 = exp(-lp2);
          
      double res, err;
      gsl_integration_qags(&f, xp2, xp1, 1.0e-8, 1.0e-8, wsize, w, &res, &err);
      ostr << l1  << " " << 0.5*(l1+l2)   << " " << l2  << " " 
           << lp1 << " " << 0.5*(lp1+lp2) << " " << lp2 << " " 
           << res/((x1-x2)*(xp1-xp2)) << " " << err/((x1-x2)*(xp1-xp2)) << endl;
    }
    ostr << endl;
  }
  
  ostr << endl << endl;


  gsl_integration_workspace_free(w);
  gsl_integration_workspace_free(wnest);

  return 0;
}


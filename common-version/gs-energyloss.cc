#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

#include "event.hh"
#include "generator.hh"
#include "meananderr.hh"

#include <vector>
#include <algorithm>

#include <CmdLine.hh>

using namespace std;

int main(int argc, char *argv[]){
  CmdLine cmd(argc, argv);

  // get the run parameters
  double epsilon = cmd.value("-epsilon", 1e-8);
  double xmin    = cmd.value("-xmin",    1e-3);
  unsigned int seed = cmd.value<unsigned int>("-seed", 1);
  unsigned int nev  = cmd.value<unsigned int>("-nev", 1000);
  string out = cmd.value<string>("-out");
  double tmax = cmd.value("-tmax", 3.0);
  double dt   = cmd.value("-dt",   0.1);

  cmd.assert_all_options_used();

  // initialise the generators

  GeneratorInMediumSimple g_simple;
  g_simple.set_seed(seed);
  g_simple.set_recursive_branching(false);

  GeneratorInMedium g_full;
  g_full.set_seed(seed);
  g_full.set_recursive_branching(false);

  // structures for the results
  unsigned int nt = (unsigned int)(tmax/dt+0.1);

  vector<MeanAndErr> energy_simple(nt);
  vector<MeanAndErr> energy_full(nt);
  vector<double> expected_simple_mean(nt);
  vector<double> expected_simple_variance(nt);

  // computed expected analytic results for the simple kernel
  for (unsigned int it=0;it<nt;++it){
    double tau = (it+1)*dt;
    expected_simple_mean[it] = exp(-M_PI*tau*tau);
    double erf1=gsl_sf_erf (sqrt(M_PI)*tau);
    double erf2=gsl_sf_erf (2*sqrt(M_PI)*tau);
    expected_simple_variance[it] = 2*M_PI*tau*(erf1-erf2)+2*exp(-M_PI*tau*tau)
      - exp(-4*M_PI*tau*tau) - exp(-2*M_PI*tau*tau);
  }

  // event loop
  unsigned int nev_save=100;
  for (unsigned int iev = 0; iev < nev; ++iev ){
    // generate th esuimple and full events
    g_simple.generate_event(tmax,epsilon,xmin);
    g_full  .generate_event(tmax,epsilon,xmin);

    // compute the remaining energy for the simple event    
    vector<double> xsums(nt, 0.0);
    for (const Particle & p : g_simple.event().particles()){ 
      if (p.x()<xmin) continue;
      double t0 = p.start_time();
      double t1 = p.end_time();

      unsigned int itmin = (unsigned int)(t0/tmax * nt);
      unsigned int itmax = (t1<0) ? nt : int(t1/tmax * nt);
      assert(itmax <= nt);

      for (unsigned int it=itmin; it<itmax; ++it) xsums[it] += p.x();
    }
    for (unsigned int it=0; it<nt; ++it)
      energy_simple[it].addEntry(xsums[it]);

    // now do the same for the full event
    for (unsigned int it=0; it<nt; ++it) xsums[it] = 0.0;
    for (const Particle & p : g_full.event().particles()){ 
      if (p.x()<xmin) continue;
      double t0 = p.start_time();
      double t1 = p.end_time();

      unsigned int itmin = (unsigned int)(t0/tmax * nt);
      unsigned int itmax = (t1<0) ? nt : int(t1/tmax * nt);
      assert(itmax <= nt);

      for (unsigned int it=itmin; it<itmax; ++it) xsums[it] += p.x();
    }
    for (unsigned int it=0; it<nt; ++it)
      energy_full[it].addEntry(xsums[it]);

    // periodic output
    if ((iev+1) % nev_save == 0){
      cout << "Output for nev = " << iev+1 << endl;

      ofstream ostr(out.c_str());
      ostr << "# ran: " << cmd.command_line() << endl;
      ostr << "# nev = " << iev+1 << endl;
      ostr << "# t  simple_avg_th  simple_var_th  simple_avg_num  simple_var_num  full_avg_num  full_var_num" << endl;
      for (unsigned int it=0; it<nt; ++it){
        ostr << (it+1)*dt << "  "
             << expected_simple_mean[it] << "  "
             << expected_simple_variance[it] << "  "
             << energy_simple[it].mean() << "  "
             << energy_simple[it].variance() << "  "
             << energy_full[it].mean()<< "  "
             << energy_full[it].variance() << endl;
      }
      if ((iev+1)==15*nev_save) nev_save*=10;
    } // output
  } // event loop

  ofstream ostr(out.c_str());
  ostr << "# ran: " << cmd.command_line() << endl;
  ostr << "# nev = " << nev << endl;
  ostr << "# t  simple_avg_th  simple_var_th  simple_avg_num  simple_var_num  full_avg_num  full_var_num" << endl;
  for (unsigned int it=0; it<nt; ++it){
    ostr << (it+1)*dt << "  "
         << expected_simple_mean[it] << "  "
         << expected_simple_variance[it] << "  "
         << energy_simple[it].mean() << "  "
         << energy_simple[it].variance() << "  "
         << energy_full[it].mean()<< "  "
         << energy_full[it].variance() << endl;
  }

  return 0;
}

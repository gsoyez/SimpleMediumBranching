#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "event.hh"
#include "generator.hh"
#include <CmdLine.hh>
#include <SimpleHist.hh>
#include "SimpleHist2D.hh"

using namespace std;

  

//========================================================================
// save results to file
void write(const string &out, vector<SimpleHist> &distribs, 
           SimpleHist &D2, SimpleHist2D &D2full,
           ostringstream &header, unsigned int nev, double tmax){
  ofstream ostr(out.c_str());

  ostr << header.str();
  ostr << "#" << endl;
  ostr << "# nev = " << nev << endl;
  ostr << "#" << endl;
  
  // output the results
  ostr << "# dN/dlogx" << endl;
  ostr << "#columns are log(1/x) dN/dx  error" << endl;

  // D's for slices in time
  for (unsigned int it=0; it<distribs.size(); it++){
    ostr << "# D(x,t=" << (it+1.0)/distribs.size()*tmax << endl;
    // construct the error
    SimpleHist err = sqrt(distribs[it]);
    output(distribs[it], err, &ostr, 1.0/(distribs[it].binsize()*nev));
    ostr << endl << endl;
  }

  // D2 (diagonal part) at fimal time
  ostr << "# D2 (final time)" << endl;
  SimpleHist err = sqrt(D2);
  output(D2, err, &ostr, 1.0/(D2.binsize()*nev));
  ostr << endl << endl;

  // D2, full 2D scane, final time
  ostr << "# D2full (final time)" << endl;
  SimpleHist2D err2 = sqrt(D2full);
  output(D2full, err2, &ostr, 1.0/(D2full.u_binsize()*D2full.v_binsize()*nev));
  ostr << endl << endl;

}

//========================================================================
// go around the fact that both generators are not derived from the
// same base class, by templating the analysis
template<typename GeneratorType>
class Analysis{
public:
  Analysis(CmdLine &cmd){
    // parse command line
    tmax = cmd.value("-tmax", 1.0);
    eps  = cmd.value("-eps",  1.0e-9);
    xmin = cmd.value("-xmin", 1.0e-4);
    nev  = cmd.value<unsigned int>("-nev", 1000);
    rseq = cmd.value<unsigned int>("-rseq", 1);
    nbin = cmd.value<unsigned int>("-nbin", 160);
    nD2  = cmd.value<unsigned int>("-nD2", 40);
    nt   = cmd.value<unsigned int>("-nt", 10);
    out = cmd.value<string>("-out");    

    //----------------------------------------------------------------------
    // a header with or setup
    header << "# Ran: " << cmd.command_line() << endl;
    header << "#" << endl;
    header << "# tmax = " << tmax << endl;
    header << "# eps  = " << eps  << endl;
    header << "# xmin = " << xmin << endl;
    header << "# rseq = " << rseq << endl;
  }

  void run(){

    //----------------------------------------------------------------------
    // setup things needed for the event generation
    GeneratorType gen(rseq);
    vector<SimpleHist> distribs(nt, SimpleHist(0.0, log(1.0/xmin), nbin));
    SimpleHist D2(0.0, log(1.0/xmin), nbin);
    SimpleHist2D D2full(0.0, log(1.0/xmin), nD2, 0.0, log(1.0/xmin), nD2);

    //----------------------------------------------------------------------
    // event loop
    unsigned int nsave = 10;
    for (unsigned int iev = 0; iev<nev; ++iev){
      gen.generate_event(tmax,eps,xmin);
      
      // we only record the final particles
      vector<double> x_above_xmin;
      for (const Particle & p : gen.event().particles()){ 
        if (p.x()<xmin) continue;
        if (p.is_final())
          x_above_xmin.push_back(p.x());

        double t0 = p.start_time();
        double t1 = p.end_time();

        unsigned int itmin = int(t0/tmax * nt);
        unsigned int itmax = (t1<0) ? nt : int(t1/tmax * nt);
        assert(itmax <= nt);

        for (unsigned int it=itmin; it<itmax; ++it)
          distribs[it].add_entry(log(1.0/p.x()));
      }

      // now bin D2 
      for (unsigned int i1 = 0; i1+1 < x_above_xmin.size(); ++i1){
        const double &lx1 = log(1.0/x_above_xmin[i1]);
        unsigned int ibin1 = D2.bin(lx1);

        for (unsigned int i2 = i1+1; i2 < x_above_xmin.size(); ++i2){
          const double &lx2 = log(1.0/x_above_xmin[i2]);
          D2full.add_entry(lx1, lx2);
          D2full.add_entry(lx2, lx1);
          if (D2.bin(lx2)==ibin1) D2.add_entry(lx1);
        }
      }
      
            
      if (((iev+1) % nsave == 0)){
        cout << "Output for nev = " << iev+1 << endl;
        write(out, distribs, D2, D2full, header, iev+1, tmax);
        if ((iev+1) == 15*nsave) nsave*=10;
      }
    }
    
    write(out, distribs, D2, D2full, header, nev, tmax);
  }
  
protected:
  double tmax;       ///< max time we evolve to
  double eps;        ///< cutoff in x for emissions
  double xmin;       ///< cutoff in x for branchings
  unsigned int nev ; ///< number of events
  unsigned int rseq; ///< randon sequence
  unsigned int nbin; ///< number of bins of x
  unsigned int nD2;  ///< number of bins for the 2D binning of D2
  unsigned int nt;   ///< number of slices in t (1...nt * (tmax/nt))
  string out;        ///< where to save the results
  ostringstream header;
};

//========================================================================
int main(int argc, char *argv[]){
  //----------------------------------------------------------------------
  // parse the command line
  CmdLine cmd(argc, argv);
  
  if (cmd.present("-simple")){
    Analysis<GeneratorInMediumSimple> analysis(cmd);
    cmd.assert_all_options_used();
    analysis.run();
    return 0;
  }

  Analysis<GeneratorInMedium> analysis(cmd);
  cmd.assert_all_options_used();
  analysis.run();
  return 0;
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "event.hh"
#include "generator.hh"
#include <CmdLine.hh>
#include <SimpleHist.hh>

using namespace std;

  

//========================================================================
// save results to file
void write(const string &out, SimpleHist &distrib, ostringstream &header, unsigned int nev){
  ofstream ostr(out.c_str());

  ostr << header.str();
  ostr << "#" << endl;
  ostr << "# nev = " << nev << endl;
  ostr << "#" << endl;
  
  // construct the error
  SimpleHist err = sqrt(distrib);

  // output the results
  ostr << "# dN/dlogx" << endl;
  ostr << "#columns are log(1/x) dN/dx  error" << endl;
  output(distrib, err, &ostr, 1.0/(distrib.binsize()*nev));
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
    SimpleHist distrib(0.0, log(1.0/eps), nbin);

    //----------------------------------------------------------------------
    // event loop
    unsigned int nsave = 10;
    for (unsigned int iev = 0; iev<nev; ++iev){
      gen.generate_event(tmax,eps,xmin);
      
      // we only record the final particles
      for (double x : gen.event().final_particles()){
        distrib.add_entry(log(1.0/x));
      }
      
      if (((iev+1) % nsave == 0)){
        cout << "Output for nev = " << iev+1 << endl;
        write(out, distrib, header, iev+1);
        if ((iev+1) == 15*nsave) nsave*=10;
      }
    }
    
    write(out, distrib, header, nev);
  }
  
protected:
  double tmax;       ///< max time we evolve to
  double eps;        ///< cutoff in x for emissions
  double xmin;       ///< cutoff in x for branchings
  unsigned int nev ; ///< number of events
  unsigned int rseq; ///< randon sequence
  unsigned int nbin; ///< number of bins of x
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

#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <cstdlib>
#include <gsl/gsl_rng.h>

using namespace std;

/// Useful constants

namespace QCD
{
    // Color coefficients
    const double TR = 1./2.;
    const double CA = 3;
    const double CF = 4./3;

    // Coupling constant
    const double alphaS = 0.1; ///< Strong coupling constant (taken independent of the energy scale)
    const double alphaSbare = alphaS/(2*M_PI); ///< Normalized strong coupling constant
}

using namespace QCD;

const double epsilon = 0.000001; /// Minimal z value for detectable outgoing gluon
const double time_max = 10; /// Minimal value for the angle of a resolvable outgoing gluon
const double time_min = 0; /// Maximal value for the angle, roughly the width of the jet



///##########################################################################################################################


///                                 Defining splitting functions


// Pgg splitting kernel as a function of z.
double Pgg(double z)
{
	return CA*(pow(z,4)+1+pow(1-z,4))/(z*(1-z));
}

// Integration of the Pgg kernel over detectable energy fraction, between epsilon/x=z_min and 1/2.
double IntPgg(double z_min)
{
    // z_min=epsilon/xparent
    return -CA*(11./6.+z_min*(-4+z_min*(1-2./3.*z_min))-4*atanh(1-2*z_min));
}

// Pqqbar splitting kernel as a function of z.
double Pqqbar(double z)
{
    return TR*(pow(z,2)+pow(1-z,2));
}

// Integration of the Pqqbar kernel over detectable energy fraction, between epsilon/x=z_min and 1/2.
double IntPqqbar(double z_min)
{
    return TR*(1./3.+z_min*(-1+(z_min*(1-2./3.*z_min))));
}

// Pqq splitting kernel as a function of z.
double Pqq(double z)
{
    return CF*(1+pow(z,2))/(1-z);
}

// Integration of the Pqqbar kernel between z_min and 1-z_min.
double IntPqq(double z_min)
{
    return CF*(-3./2.+3*z_min+4*atanh(1-2*z_min));
}

///##########################################################################################################################


///                 Random variable generator following the Sudakov probability distribution


// Inverse of the Sudakov probability density
double Inverse_Sudakov(double u, double time_parent, double xparent)
{
    return time_parent-1/(alphaSbare*IntPgg(epsilon/xparent))*log(u);
}

// Random variable generator following the Sudakov probability density
double alea_generator_sudakov(double time_parent, double xparent,gsl_rng *r)
{
	double alea;
	alea = gsl_rng_uniform(r);
	return Inverse_Sudakov(alea,time_parent,xparent);
}


///##########################################################################################################################


///             Random variable generator following the Pgg probability density


// Approximate Pgg density between 1/2 and 1-epsilon for the veto algoritm
double alea_generator_Pgg_approx(double xparent,gsl_rng *r)
{
	double norm, al;
	norm = 2*log(-1+1/(epsilon/xparent)); /// Integral of 2/(z(1-z)) between epsilon/x and 1/2
	al = gsl_rng_uniform(r)*norm;
	return 1/(1+exp(-al/2.)); /// Generates a number between 1/2 and 1-epsilon/x
}

// Implementation of the veto algorithm, Pgg density between 1/2 and 1-epsilon
double alea_generator_Pgg(double xparent,gsl_rng *r)
{
	double output, z_candidate, alea;
	while(1)
    {
		z_candidate = alea_generator_Pgg_approx(xparent,r);
		alea = gsl_rng_uniform(r)*CA*2/(z_candidate*(1-z_candidate));
		output = z_candidate;
		if(alea<=Pgg(z_candidate))
        {
			break;
        }
    }
	return output; /// Generates a number between 1/2 and 1-epsilon/
}



///##########################################################################################################################


///              Random variable for the medium induced K(z) probability density


// Approximate K density between epsilon and 1/2
double alea_generator_K_approx(double xparent, gsl_rng *r)
{
    double z_min = epsilon/xparent;
    double norm;
    norm = (2*(1-2*z_min))/sqrt(z_min*(1-z_min)); /// Integral of 1/(z*(1-z)0^(3/2) between epsilon/x and 1/2
    double al = gsl_rng_uniform(r)*norm;
    al*=al;
    return (16+al-sqrt(16*al+al*al))/(2*(16+al)); /// Generates a number between epsilon/x and 1/2
}


///##########################################################################################################################


///                   Generating one cascade


// At each step, for all the gluons produced at the precedent step,
// an angle and a energy fraction z are chosen randomly. If the angle
// is above the minimal resolution angle and if the outgoing gluons
// carry and energy above the threshold of detection, two gluons and
// a new vertex are added to the two vectors "partons" and "vertices".

// The function returns the list of all the gluons produced during the
// cascade (not only the final ones).



class Parton
{
    public:

    // Constructors
    Parton()
    {this->starting_vertex=-1;this->ending_vertex=0; this->x=1;this->time=0;};
    Parton(int starting_vertex,int ending_vertex, double x, double time)
    {this->starting_vertex=starting_vertex; this->ending_vertex=ending_vertex; this->x=x;this->time=time;};

    // Getters
    int get_starting_vertex() const { return this->starting_vertex;}
    double get_ending_vertex() const { return this->ending_vertex;}
    double get_x() const {return this-> x;}
    double get_time() const {return this->time;}

    // Setters
    void set_starting_vertex(int starting_vertex) {this->starting_vertex=starting_vertex;}
    void set_ending_vertex(int starting_vertex) {this->starting_vertex=starting_vertex;}
    void set_x(double x) {this->x=x;}
    void set_time(double time) {this->time=time;}

	private:


	int starting_vertex;   ///< Index in given list of the vertex which created the parton
	double ending_vertex;  ///< Index in given list of the vertex which deleted the parton
	double x;              ///< Total fraction of energy x carried by the parton
	double time;           ///< Useless !
};


class Vertex
{
    public:
    // Constructors
    Vertex()
    {this->index_parent=0; this->index_first=1; this->index_second=2; this->time=0;};
    Vertex(int index_parent, int index_first, int index_second, double time)
    {this->index_parent=index_parent; this->index_first=index_first; this->index_second=index_second; this->time=time;}

    // Getters
    int get_index_parent() const { return this->index_parent;}
    int get_index_first() const { return this->index_first;}
    int get_index_second() const { return this->index_second;}
    double get_time() const { return this->time;}

    // Setters
    void set_index_parent(int index_parent) {this->index_parent=index_parent;}
    void set_index_first(int index_first) {this->index_first=index_first;}
    void set_index_second(int index_second) {this->index_second=index_second;}
    void set_time(double time) {this->time=time;}

	private:

    int index_parent;  ///< Index in a given list of the parent parton.
	int index_first;   ///< Index in a given list of the first parton going out the vertex
	int index_second;  ///< Index in a given list of the second parton going out the vertex
	double time;       ///< Square of the angle associated with the vertex
};

//class Event
//{
//    public:
//    vector<Parton> parton_produced;
//    vector<Vertex> vertex_produced;
//};


vector<Parton>  PartonShower(double time_ini, double x_ini, gsl_rng *r)
{

	// Initializing the first parton
	vector<Parton> all_partons(1);
	Parton first_parton(-1,0,x_ini,0);
	all_partons[0]=first_parton;
	vector<Vertex> all_vertices(0);

    int i=0;
    int N_vertex = -1;

	while(i<all_partons.size())
    {
        // Angle and total energy fraction x carried by the precedent gluon
        double last_time,last_x;
        int last_vertex;
        last_time = all_partons[i].get_time();
        last_x = all_partons[i].get_x();
        last_vertex = all_partons[i].get_ending_vertex();

        if(last_time<time_max && last_x>2*epsilon)
        {
            N_vertex++;

            double new_time, new_z;
            new_time = alea_generator_sudakov(last_time,last_x,r);
            new_z = alea_generator_Pgg(last_x,r);

            Parton new_parton1(last_vertex,N_vertex,new_z*last_x,new_time);
            Parton new_parton2(last_vertex,N_vertex,(1-new_z)*last_x,new_time);
            all_partons.push_back(new_parton1);
            all_partons.push_back(new_parton2);

            Vertex new_vertex(i,all_partons.size()-2,all_partons.size()-1,new_time);
            all_vertices.push_back(new_vertex);
        }
        i++;
    }
    return all_partons;
}


///##########################################################################################################################

// The main program provides an output file with the total
// number of gluons at the final step in each (logarithmic)
// bin.


int main()
{
    // Initializing the uniform random generator
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

	int N_segment; // Number of bins
	N_segment = 100;

	int N_montecarlo; // Number of cascades generated
	N_montecarlo = 10000;

	vector<int> histogramme(N_segment,0); // List with the number of gluons in each bin (to plot an histogram diagram)
	int i;
	for(i=0;i<N_montecarlo;i++)
    {
		vector<Parton> one_cascade; // Generating one cascade
		one_cascade = PartonShower(time_min,1.,r);
		int N_one_cascade;
		N_one_cascade = one_cascade.size();

		int j;
		for(j=0;j<N_one_cascade;j++)
        {
			if(one_cascade[j].get_time()>time_max && one_cascade[j].get_x()>=epsilon && one_cascade[j].get_x() <= 1-epsilon) // Selects the gluons produced at the final step of the shower
            {
				double final_x;
				int k_segment;
				final_x = one_cascade[j].get_x();
				k_segment = floor(-log(final_x/epsilon)*N_segment/log(epsilon)); // Finds the bin (with logarithmic scale)
				histogramme[k_segment] =histogramme[k_segment]+1; // Adds one gluon to this bin
            }
        }
    }


	// The vector "histogramme" is given in a file "output.txt"
	fstream histo;
	histo.open("output.dat",ios::out);
	for(i=0;i<N_segment;i++)
    {
		histo << epsilon*exp(-i*log(epsilon)/N_segment);
		histo << " ";
		histo << histogramme[i]*1./N_montecarlo*1./(-log(epsilon)/N_segment) << endl;
    }
	histo.close();
	system("gnuplot plot.gnu");
	return 0;
}

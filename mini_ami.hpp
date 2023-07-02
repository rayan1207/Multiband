#ifndef MINI_AMI_HPP
#define MINI_AMI_HPP


#include "ami_base.hpp"
#include "ami_calc.hpp"
#include "amigraph.hpp"
#include <cassert>
#include <experimental/filesystem>
#include <cmath>
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <unistd.h>
#include <mpi.h>
#include <chrono>
#include <thread>
#include <sstream>

#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <typeinfo> 
#include <string>
#include <utility> 
#include <locale>
#include <stdexcept>

class mband{
private:

 
 
public: 
std::vector<vector<int>> interaction_legs; //= {{1,1,1,1},{2,2,1,1},{2,1,2,1},{1,2,2,1},{2,1,1,2},{1,2,1,2},
//{1,1,2,2},{2,2,2,2}};

std::vector<double> int_values;//={0.69907,0.664234,0.1815454,0.1815454,0.1815454,0.1815454,0.664234,0.674536};
std::vector<double> energy;
bool Hartee_fock;
mband(std::vector<std::vector<int>> _interaction_legs,  std::vector<double> _int_values,std::vector<double> _energy, bool _Hartee_fock );

struct sampler_collector {
    AmiGraph::graph_t graph;
	AmiGraph::edge_vector_t fermionic_edge;
	std::vector<std::vector<int>> fermionic_edge_species;
	std::vector<std::vector<std::vector<int>>> interaction_species ;
	std::vector<std::vector<int>> external_line;
	std::vector<AmiBase::epsilon_t> Epsilon;
	std::vector<AmiBase::alpha_t> Alpha;
	std::vector<AmiBase::alpha_t> bosonic_Alpha;
	std::vector<std::vector<int>>  Uindex;
	
};
struct output_collector {
    std::vector<std::complex < double >> result_vec;
    std::vector<double> beta_vec;
    std::vector<double> mfreq_vec;
    std::vector<std::vector<int>> extline_vec;
    std::vector<std::vector<int>> Uindex_vec;
};
struct params_param {
	int  min_ord;		
	int  max_ord;
    int molecular;
    double E_reg;
    int lattice_type;
	int molecular_type;
    double set_precision;
    bool hatree_fock;
    double tp;
    double tperp;
	double tperp_p;
	int MC_num;
	int in;
	int out;
	double molec_beta;
	int molec_mfreq;
	int time;
	double V;
	
};


double Bilayer_Hubbard_Energy(NewAmiCalc::ext_vars ext,std::vector<double> momenta, int species,mband::params_param param);
void assign_label(AmiGraph::graph_t &g1, AmiGraph::edge_vector_t edge,vector<int> vector);
std::vector<std::vector<int>>  findmatch(std::vector<int> v1,AmiGraph::edge_vector_t &inter_vec);
std::vector<int> generate_edge_species(AmiGraph::graph_t &g, AmiGraph::edge_vector_t edge);
std::vector<int> external_species(AmiGraph::graph_t &graph);
void solve_multiband_5(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs);
void solve_multiband_4(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs);
void solve_multiband_3(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs);
void solve_multiband_2(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs);
void solve_multiband_1(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs);
//void solve_multiband_234(AmiGraph::graph_t graph,std::vector<std::vector<int>> interaction_leg,int ord);
void print_match(std::vector<int> vec,int ord);
void reset_species(AmiGraph::graph_t &graph,std::vector<AmiGraph::edge_vector_t> int_vector);
void print_assigned_species(std::vector<std::vector<std::vector<int>>> all_species);
void find_interaction(AmiGraph::graph_t &graph, AmiGraph::edge_vector_t &b_vector, std::vector<AmiGraph::edge_vector_t> &f_vector);
void print_interactions(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t b_vector, std::vector<AmiGraph::edge_vector_t> f_vector);
void  generate_eps_alpha(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge, std::vector<std::vector<int>> &ept, 
std::vector<std::vector<int>>  &alpha);


//std::vector<int>  interaction_index(std::vector<std::vector<int>> int_species);
std::vector<int> interaction_index(const std::vector<std::vector<int>>& int_species);

//template<typename T>
//void print2d(const std::vector<std::vector<T>>& vec);
void solve_multiband(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs);
void molecular_solver_ext( AmiGraph::graph_t &gself, output_collector& collector,std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec, std::vector<int> line,std::string outputfilem );
void molecular_solver( AmiGraph::graph_t &gself, output_collector& collector,std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec,std::string outputfile );
void sigma_sampler( AmiGraph::graph_t &gself, sampler_collector& collector);
void write_output(std::string outputfile,output_collector& collector,std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec);
void calculate_sampled_sigma(AmiGraph::graph_t &gself, sampler_collector& samp_collector,  output_collector& out_collector, std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec );
void calculate_sampled_sigma_ext(AmiGraph::graph_t &gself, sampler_collector& samp_collector,  output_collector& out_collector, std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec,std::vector<int> line );
std::pair<std::complex<double>, std::complex<double>> lcalc_sampled_sigma(AmiGraph::graph_t &gself, std::vector<AmiBase::epsilon_t>& Epsilon, std::vector<AmiBase::alpha_t>& Alpha,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<int> &Utype, std::vector<int>& Species,NewAmiCalc::ext_vars& ext_params,int MC_num,params_param& param);


double non_local_U_formfactor(std::vector<double> vq);
std::vector<std::vector<int>> find_non_local_bosonic_alpha(std::vector<std::vector<int>> bosonic_Alpha, std::vector<int> Uindex);

std::vector<double> sumVectors(std::vector<std::vector<double>> vectors);
std::vector<std::complex<double>> convertToComplex(const std::vector<double> vec);
std::vector<std::vector<double>>  band_to_hab(std::vector<std::vector<int>> band);
std::vector<std::complex<double>> generate_ept(std::vector<std::vector<int>> epsilon, std::vector<double> band_value);
std::vector<AmiBase::epsilon_t> updateEpsilon(const std::vector<AmiBase::epsilon_t>& epsilon, const std::vector<double>& energy);
double perturb(double number);
double Stdev(double total_sq, double mean, int n);
//double Umatch(std::vector<std::vector<int>> int_matrix, std::vector<double> int_value, std::vector<std::vector<int>> int_species);
double Umatch(const std::vector<std::vector<int>>& int_matrix, const std::vector<double>& int_value, 
                     const std::vector<std::vector<int>>& int_species);

std::vector<int> Hartee_fock_filter(AmiGraph::edge_vector_t &fermionic_edge);
void filter(std::vector<std::vector<int>>& possible_species, const std::vector<int>& list);

//////////////////////////////////lattice_stuff///////////////////////////

double Hubbard_Energy(NewAmiCalc::ext_vars ext,std::vector<double> momenta, int species,params_param param);
//std::complex<double> lcalc_sampled_sigma(AmiGraph::graph_t &gself, std::vector<AmiBase::epsilon_t>& Epsilon, std::vector<AmiBase::alpha_t>& Alpha, std::vector<int>& Species,
//    NewAmiCalc::ext_vars& ext_params, int MC_num, int lattice_type);
};
/// A few simple functions
template<typename T>
inline void print2d( std::vector< std::vector<T>> vec)
{
    for ( auto row : vec) {
        for ( auto elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
	
}


template<typename T> 
inline void print1d( std::vector<T>& vec) {
  std::cout << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    std::cout << vec[i];
    if (i != vec.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << "]\n";
}




inline std::vector<std::vector<int>> readFile(const char* filename) {
    std::ifstream inFile(filename);
    std::string line;
    std::vector<std::vector<int>> data;

    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        std::vector<int> row;
        int val;

        for (int i = 0; i < 4; i++) {
            iss >> val;
            row.push_back(val);
        }

        data.push_back(row);
    }

    return data;
}


inline std::vector<double> readFile1(const char* filename, int n) {
    std::ifstream inFile(filename);
    std::string line;
    std::vector<double> columnValues;

    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        double val;

        for (int i = 1; i < n; i++) {
            if (iss >> val) {
               
            } else {
                std::cerr << "Error: Could not read column " << n << " from file " << filename << std::endl;
                return std::vector<double>();
            }
        }

        if (iss >> val) {
            
            columnValues.push_back(val);
        } else {
            std::cerr << "Error: Could not read column " << n << " from file " << filename << std::endl;
            return std::vector<double>();
        }
    }

    return columnValues;
}
std::string trim(const std::string& str);
void parseLine(const std::string& line, std::string& paramName, std::string& paramValue);
void params_loader(const std::string& filename, mband::params_param& params);
#endif // MINI_AMI_HPP


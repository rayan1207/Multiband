
#include "ami_base.hpp"
#include "ami_calc.hpp"
#include "amigraph.hpp"
#include <cassert>
#include <experimental/filesystem>

#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <unistd.h>

#include <chrono>
#include <thread>

#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>


class mband{
private:

 
 
public: 
std::vector<vector<int>> interaction_legs; //= {{1,1,1,1},{2,2,1,1},{2,1,2,1},{1,2,2,1},{2,1,1,2},{1,2,1,2},
//{1,1,2,2},{2,2,2,2}};

std::vector<double> int_values;//={0.69907,0.664234,0.1815454,0.1815454,0.1815454,0.1815454,0.664234,0.674536};
std::vector<double> energy;
bool Hartee_fock;
mband(std::vector<std::vector<int>> _interaction_legs,  std::vector<double> _int_values,std::vector<double> _energy, bool _Hartee_fock );





void assign_label(AmiGraph::graph_t &g1, AmiGraph::edge_vector_t edge,vector<int> vector);
std::vector<std::vector<int>>  findmatch(std::vector<int> v1,AmiGraph::edge_vector_t &inter_vec);
std::vector<int> generate_edge_species(AmiGraph::graph_t &g, AmiGraph::edge_vector_t edge);
std::vector<int> external_species(AmiGraph::graph_t &graph);
void solve_multiband_5(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &ext_legs);
void solve_multiband_4(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &ext_legs);
void solve_multiband_3(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &ext_legs);
void solve_multiband_2(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &ext_legs);
void solve_multiband_1(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &ext_legs);
//void solve_multiband_234(AmiGraph::graph_t graph,std::vector<std::vector<int>> interaction_leg,int ord);
void print_match(std::vector<int> vec,int ord);
void reset_species(AmiGraph::graph_t &graph,std::vector<AmiGraph::edge_vector_t> int_vector);
void print_assigned_species(std::vector<std::vector<std::vector<int>>> all_species);
void find_interaction(AmiGraph::graph_t &graph, AmiGraph::edge_vector_t &b_vector, std::vector<AmiGraph::edge_vector_t> &f_vector);
void print_interactions(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t b_vector, std::vector<AmiGraph::edge_vector_t> f_vector);
void  generate_eps_alpha(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge, std::vector<std::vector<int>> &ept, 
std::vector<std::vector<int>>  &alpha);
std::vector<int>  interaction_index(std::vector<std::vector<int>> int_species);
//template<typename T>
//void print2d(const std::vector<std::vector<T>>& vec);
void solve_multiband(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &ext_legs);

///


std::vector<double> sumVectors(std::vector<std::vector<double>> vectors);
std::vector<std::complex<double>> convertToComplex(const std::vector<double> vec);
std::vector<std::vector<double>>  band_to_hab(std::vector<std::vector<int>> band);
std::vector<std::complex<double>> generate_ept(std::vector<std::vector<int>> epsilon, std::vector<double> band_value);
double perturb(double number);
double Umatch(std::vector<std::vector<int>> int_matrix, std::vector<double> int_value, std::vector<std::vector<int>> int_species);
std::vector<int> Hartee_fock_filter(AmiGraph::edge_vector_t &fermionic_edge);
void filter(std::vector<std::vector<int>>& possible_species, const std::vector<int>& list);
};



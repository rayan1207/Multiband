
#include "mini_ami.hpp"
#include <typeinfo> 
#include <vector>
#include <utility> 

void print2d( std::vector<std::vector<int>> vec)
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
void print1d( std::vector<T>& vec) {
  std::cout << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    std::cout << vec[i];
    if (i != vec.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << "]\n";
}




std::vector<std::vector<int>> readFile(const char* filename) {
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


std::vector<double> readFile1(const char* filename, int n) {
    std::ifstream inFile(filename);
    std::string line;
    std::vector<double> columnValues;

    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        double val;

        for (int i = 1; i < n; i++) {
            if (iss >> val) {
                // Read and discard values until we get to the nth column
            } else {
                std::cerr << "Error: Could not read column " << n << " from file " << filename << std::endl;
                return std::vector<double>();
            }
        }

        if (iss >> val) {
            // Successfully read the nth column value, store it in the vector
            columnValues.push_back(val);
        } else {
            std::cerr << "Error: Could not read column " << n << " from file " << filename << std::endl;
            return std::vector<double>();
        }
    }

    return columnValues;
}


int main(int argc, char** argv)
{ 


    AmiBase ami;
	std::vector<vector<int>> xx = {{1,1,1,1},{2,2,1,1},{2,1,2,1},{1,2,2,1},{2,1,1,2},{1,2,1,2},{1,1,2,2},{2,2,2,2}};
    std::vector<vector<int>> hubbard = {{2,2,1,1},{1,1,2,2}};

	 //std::vector<std::vector<int>> xx = readFile("data.txt");
	 //std::vector<double> yy = readFile1("data.txt",5);
	 //std::vector<double> zz = readFile1("data_h.txt",2);
	 
	 
	 
    std::vector<double> yy={0.69907,0.664234,0.1815454,0.1815454,0.1815454,0.1815454,0.664234,0.674536};
	
	 std::vector<double> zz = {-0.5825365736653964,0.6670627411988165};
	
	

	
	std::cout<<"This will be a minimal working tutorial/example of AMI+libamigraph"<<std::endl;
	print1d(yy);
	std::this_thread::sleep_for(std::chrono::milliseconds(1800)); // Just pause for 2 seconds so user can follow 
	
	std::cout<<std::endl;
	
	
	
	//
	
	std::cout<<"First need to generate an instance of the AmiGraph library"<<std::endl;
	std::this_thread::sleep_for(std::chrono::milliseconds(150));
		int seed=3;	
		std::cout<< "Constructing AmiGraph object using seed: "<<seed <<" " <<std::endl;
		AmiGraph g(AmiBase::Sigma, seed);	

	std::string infile("ext_vars.dat");
	
		std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
		NewAmiCalc::external_variable_list extern_list;

	// read column data into a file - you can get external variables however you would like	
	g.ami.read_external(infile, extern_list);	
	
	std::cout<<"External parameters read are"<<std::endl;
for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_freq_[0]<<std::endl;
		}	

// Next lets load the graphs 

	AmiGraph::gg_matrix_t ggm; 
	// ggm stands for 'graph group matrix' it is a multidimensional array separating graphs into orders and group
       
	std::cout<<"To get started, we need to load some graphs"<<std::endl;
	std::cout<<"Attempting to load self-energy graphs from ../../example_graphs"<<std::endl;
	int max=4;
	g.read_ggmp("../example_graphs/ggm_all/",ggm, max);
	std::cout<<"Completed read"<<std::endl;

	//g.mpi_print_ggm(ggm, 0);
	
	//g.ggm_remove_pairs(ggm,0); // remove the pair that was defined at third order - we don't need it for now 
	//std::cout<<"Removed the cancelling pair - now they are in a group"<<std::endl;
	g.mpi_print_ggm(ggm,0);
	
	
	std::cout<<std::endl;
	
	std::this_thread::sleep_for(std::chrono::milliseconds(1500)); 
	std::cout<<"Now lets label only the second order graphs "<<std::endl;
	
	g.ggm_label(ggm,0);
	std::cout<<"All done! That was easy!"<<std::endl;
	std::this_thread::sleep_for(std::chrono::milliseconds(1500));
	
///Rayan function	
int  ord = 2;
	AmiGraph::graph_t  gself2 = ggm[ord][0].graph_vec[0];
    //AmiGraph::graph_t  gself = ggm[ord][1].graph_vec[0];

	
	
	std::cout<<"Lets check that a graph is labelled"<<std::endl;
	g.print_all_edge_info(gself2);
    
	
	mband mb(xx,yy,zz);
	std::cout<< mb.interaction_legs.size();
	//mb.solve_multiband_2(ggm[2][0].graph_vec[0]);
	//AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species
	AmiGraph::edge_vector_t o2fedge;
	/// contains possible band index for o2fedge( internal fermionic edge) in same order as o2 fedge
    std::vector<std::vector<int>> o2fedge_species;
	
	std::vector<std::vector<std::vector<int>>> o2int_species;
	std::vector<AmiBase::epsilon_t> Epsilon;
	std::vector<AmiBase::alpha_t> Alpha;
	/// provide a second order graph in gself2, and it fills in the internal fermionic edge pointed in o2fedge,
	//all the possible species o2fedge_species(8)  in the same order o2 fedge // finally o2int_species contains interactions
	//that give rise to o2fedge-species.
	
	mb.solve_multiband(gself2,o2fedge,o2fedge_species,o2int_species,ord);
	
	

	
     	
	
	//std::cout<< " total number of cycles " << num << std::endl;
    std::cout <<"printing the band index for each internal fermionic edge" <<std::endl;
	
     for (int i = 0; i < o2fedge_species.size();i++){
		 std::cout<< "(";
		 for (int j = 0; j < o2fedge_species[i].size();j++){
			 std::cout<<o2fedge_species[i][j];	 
		}
	 std::cout<<")\n\n";

	 };
	 
	 ///filling in vectors of epsilon and alpha based on labels
	mb.generate_eps_alpha(gself2,o2fedge,Epsilon,Alpha);
	
	std::vector<std::vector<double>> kk = mb.band_to_hab(o2fedge_species);
	std::vector<std::vector<std::complex<double> >> energy_V;
	
	for (auto vec: kk){
		energy_V.push_back(mb.generate_ept(Epsilon, vec));
		//convering energies to complex format for energy_t
		//energy_V.push_back(mb.convertToComplex(vec));
	}
	 
	
      
	  
	 
	 


	 g.print_all_edge_info(gself2);	
	 
	 
	std::cout<<"Lets find the each edge and possible band indexes! \n";

	for (int i=0;i <o2fedge.size();i++){	
		g.print_edge_info(o2fedge[i],gself2);
		//ept = gself2[o2fedge[i]].g_struct_.eps_;
		std::cout<<"epsilon we generate for AMI is ";
		print1d(Epsilon[i]);
		std::cout<<" alpha is " ;
		print1d(Alpha[i]);
		std::cout<<" all possible band indexes are: ";
		for (int j = 0; j < o2fedge_species.size(); j++){
			std::cout<< o2fedge_species[j][i] <<" ";
		}
	std::cout<<std::endl <<"\n";	
	}
	std::cout << "energy V are" <<std::endl;
	for (int i = 0; i <energy_V.size();i++){
		for (int j = 0; j <energy_V[i].size();j++)
		{
			std::cout<<energy_V[i][j];
		}
		std::cout<<std::endl;
	}
	///evaluation using AMI _-_-_-_-_-*/


int n = 2*ord-1; // arbitrary number of elements in the array
AmiBase::g_struct gs[n];
std::vector<AmiBase::g_struct> gs_vec;
for(int i = 0; i < n; i++) {
     
    gs[i] = {Epsilon[i], Alpha[i]};
	gs_vec.push_back(gs[i]);
}

//AmiBase::g_struct gs[3] = {{Epsilon[0],Alpha[0]},{Epsilon[1],Alpha[1]},{Epsilon[2],Alpha[2]} };


AmiBase::g_prod_t R0 =gs_vec;
//print2d({gs[0].eps_,gs[0].alpha_});
//print2d({gs[1].eps_,gs[1].alpha_}); 
//print2d({gs[2].eps_,gs[2].alpha_}); 
AmiBase::S_t S_array;
AmiBase::P_t P_array;
AmiBase::R_t R_array;
double beta = 50;
double E_REG=0.0000000; 
int N_INT=ord;
  
AmiBase::ami_parms test_amiparms(N_INT, E_REG);

ami.construct(test_amiparms , R0 , R_array , P_array , S_array ); 

AmiBase:: frequency_t frequency;
for (int i= 0; i<ord;i++){
frequency.push_back(std::complex<double>(0,0));}
frequency.push_back(std::complex<double>(0,M_PI/beta));
std::complex<double> final_result = {0,0};
for (auto energy: energy_V){
int i =0;
//std::cout<<"Using energy_t" <<std::endl;	
for (auto e : energy){std::cout << e ;}

AmiBase::ami_vars external (energy,frequency,beta  );
std::complex < double > calc_result = ami.evaluate(test_amiparms,R_array ,P_array,S_array,external);
std::cout<<"result is " <<  calc_result*mb.Umatch(xx,yy,o2int_species[i]) <<"\n";
final_result = final_result+  calc_result*mb.Umatch(xx,yy,o2int_species[i]);

i++;
}
std::cout<<"Final Result is result is " << final_result <<"\n";
std::cout<<"Number of internal fermionic line is:" <<o2fedge.size()<< std::endl;
std::cout<<"different number of possible species arrangements are: " <<o2fedge_species.size() <<std::endl;
std::cout<<"total number possible U_abcd interaction printed above: " <<o2int_species.size() <<"\n \n";





	
	
	

		
	
	
}

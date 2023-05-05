
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


int main(int argc, char** argv)
{ 


    AmiBase ami;
	
	
     std::vector<std::vector<int>> interaction = readFile("data.txt");
	 std::vector<double> interaction_value = readFile1("data.txt",5);
	 std::vector<double> band_energy = readFile1("data_h.txt",2);
	
	 /*
	 std::vector<vector<int>> interaction = {{1,1,1,1},{2,2,1,1},{2,1,2,1},{1,2,2,1},{2,1,1,2},{1,2,1,2},{1,1,2,2},{2,2,2,2}};
     std::vector<vector<int>> hubbard = {{2,2,1,1},{1,1,2,2}};
     std::vector<double> interaction_value={0.69907,0.664234,0.1815454,0.1815454,0.1815454,0.1815454,0.664234,0.674536};
	 //std::vector<double> yy={0.69907,0.4825,0.1815,-0.4825,-0.4825,0.1815,0.4825,0.6745};
	 std::vector<double> band_energy = {0.5825365736653964,-0.6670627411988165};
	*/

		int seed=2;	
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

	//g.mpi_print_ggm(ggm,0);
	
	
	std::cout<<std::endl;
	
	std::this_thread::sleep_for(std::chrono::milliseconds(1500)); 
	std::cout<<"Now lets label only the second order graphs "<<std::endl;
	
	g.ggm_label(ggm,0);
	
	std::cout<<"All done! That was easy!"<<std::endl;
	std::this_thread::sleep_for(std::chrono::milliseconds(1500));
	
///Rayan function	
int  ord = 3;
	AmiGraph::graph_t  gself = ggm[ord][0].graph_vec[0];
    //AmiGraph::graph_t  gself = ggm[ord][1].graph_vec[0];
    
	
	
	std::cout<<"Lets check that a graph is labelled"<<std::endl;
	g.print_all_edge_info(gself);
    
	bool hf = true;
	//construction 
	mband mb(interaction,interaction_value,band_energy,hf);
	std::cout<< mb.interaction_legs.size();
	
	AmiGraph::edge_vector_t fermionic_edge;
	/// contains possible band index for fermionic_edge( internal fermionic edge) in same order as o2 fedge
    std::vector<std::vector<int>> fermionic_edge_species;
	std::vector<std::vector<std::vector<int>>> interaction_species ;
	
	std::vector<std::vector<int>> external_line;
	/// provide a second order graph in gself, and it fills in the internal fermionic edge in fermionic_edge vector,
	//all the possible species fermionic_edge_species corresponding  to fermionic_edge vector
	// finally interaction_species  contains interactions that give rise to fermionic_edge-species.
	//
	
	mb.solve_multiband(gself,fermionic_edge,fermionic_edge_species,interaction_species ,external_line);
	std::this_thread::sleep_for(std::chrono::milliseconds(1500)); 
    std::cout <<"printing the band index for each internal fermionic edge" <<std::endl;
	
     for (int i = 0; i < fermionic_edge_species.size();i++){
		 std::cout<< "(";
		 for (int j = 0; j < fermionic_edge_species[i].size();j++){
			 std::cout<<fermionic_edge_species[i][j];	 
		}
	 std::cout<<")\n\n";

	 };
	 
	 ///filling in vectors of epsilon and alpha based on labels
	std::vector<AmiBase::epsilon_t> Epsilon;
	std::vector<AmiBase::alpha_t> Alpha;
	mb.generate_eps_alpha(gself,fermionic_edge,Epsilon,Alpha);
	
	std::vector<std::vector<double>> speciesToEnergy = mb.band_to_hab(fermionic_edge_species);
	std::vector<std::vector<std::complex<double> >> energy_vector;
	//  using the epsilon labels, we put the energies in correct place. e.g. energy_vector = e1[0,1,0] + e2[1,0,0] + e3[0,0,1] where 
    // e1 may correspond to energy for band index 1, e2 energy for band index 2 and so on/	
	for (auto vec: speciesToEnergy){
		energy_vector.push_back(mb.generate_ept(Epsilon, vec));

	}
	 g.print_all_edge_info(gself);	
	 
	 
	std::cout<<"Lets find the each edge and possible band indexes! \n";

	for (int i=0;i <fermionic_edge.size();i++){	
		g.print_edge_info(fermionic_edge[i],gself);
		std::cout<<"epsilon we generate for AMI is ";
		print1d(Epsilon[i]);
		std::cout<<" alpha is " ;
		print1d(Alpha[i]);
		std::cout<<" all possible band indexes are: ";
		for (int j = 0; j < fermionic_edge_species.size(); j++){
			std::cout<< fermionic_edge_species[j][i] <<" ";
		}
	std::cout<<std::endl <<"\n";	
	}
	std::cout << "energy V are" <<std::endl;
	for (int i = 0; i <energy_vector.size();i++){
		for (int j = 0; j <energy_vector[i].size();j++)
		{
			std::cout<<energy_vector[i][j];
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
double E_REG=1e-7; 
int N_INT=ord;
  
AmiBase::ami_parms test_amiparms(N_INT, E_REG);
//ami.drop_matsubara_poles = false;
ami.construct(test_amiparms , R0 , R_array , P_array , S_array ); 

AmiBase:: frequency_t frequency;
for (int i= 0; i<ord;i++){
frequency.push_back(std::complex<double>(0,0));}
frequency.push_back(std::complex<double>(0,M_PI/beta));
std::complex<double> final_result = {0,0};
//Self energy for all different band indexes are evaluated for first matsubara frequency. This is meant as a test.
for (int i = 0; i<energy_vector.size();i++){	

AmiBase::ami_vars external (energy_vector[i],frequency,beta  );
std::complex < double > calc_result = ami.evaluate_otf(test_amiparms,R_array ,P_array,S_array,external);
std::cout<<"(" << i+1 <<")"<<  calc_result*mb.Umatch(interaction,interaction_value,interaction_species [i]) <<std::endl;
//std::cout<<"result is " <<  calc_result <<std::endl;
final_result = final_result+  calc_result*mb.Umatch(interaction,interaction_value,interaction_species [i]);
}
std::cout<<"Final Result is result is " << final_result <<"\n";
std::cout<<"Number of internal fermionic line is:" <<fermionic_edge.size()<< std::endl;
std::cout<<"different number of possible species arrangements are: " <<fermionic_edge_species.size() <<std::endl;
std::cout<<"total number possible U_abcd interaction printed above: " <<interaction_species.size() <<"\n \n";
print2d(external_line);

frequency.clear();
std::ofstream outFile("g2_band10_mfreq_beta50.txt");


for  (int m = 0; m <fermionic_edge_species.size();m++){
	for (int n=0; n < 10; n++){	
		for (int i= 0; i<ord;i++){
			frequency.push_back(std::complex<double>(0,0));}
			frequency.push_back(std::complex<double>(0,(2*n+1)*M_PI/beta));

			AmiBase::ami_vars external (energy_vector[m],frequency,beta  );
			std::complex < double > calc_result = ami.evaluate(test_amiparms,R_array ,P_array,S_array,external);//*mb.Umatch(xx,yy,interaction_species[2]) ;
			//std::cout<<"result is \n" <<  calc_result ;
			outFile << (2*n+1)*M_PI/beta << " " << calc_result.real() << " "<< calc_result.imag()<<" "<<external_line[m][0]<<" "<< external_line[m][1] <<" ";
			for (auto i : mb.interaction_index(interaction_species[m])){
				outFile<< i <<" ";				
			}
			outFile << std::endl;

	
	
frequency.clear();	
	}
}	
outFile.close();
}


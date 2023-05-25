
#include "mini_ami.hpp"

int main(int argc, char** argv)
{ 


    AmiBase ami;
	
    /*
     std::vector<std::vector<int>> interaction = readFile("three_orb_U.txt");
	 std::vector<double> interaction_value = readFile1("three_orb_U.txt",5);
	 std::vector<double> band_energy = readFile1("three_orb_E.txt",2);
	 
	 std::cout<<interaction.size();

	
	print1d(band_energy);
    */
	
    double e1 = -0.5825365736653964;
	double e2 = 0.6670627411988165;
	
	 double e1_ays  = 0.5*e1 - 0.5*e2;
	 double e2_ays  = 0.5*e2 -0.5*e1;

	std::cout <<"energies are" << e1_ays <<" " <<e2_ays;
	 std::vector<vector<int>> interaction = {{1,1,1,1},{2,2,1,1},{2,1,2,1},{1,2,2,1},{2,1,1,2},{1,2,1,2},{1,1,2,2},{2,2,2,2}};
     std::vector<vector<int>> hubbard = {{2,2,1,1},{1,1,2,2}};
     std::vector<double> interaction_value={0.69907,0.664234,0.1815454,0.1815454,0.1815454,0.1815454,0.664234,0.674536};
	 //std::vector<double> yy={0.69907,0.4825,0.1815,-0.4825,-0.4825,0.1815,0.4825,0.6745};
	 //std::vector<double> band_energy = {-0.5825365736653964,0.6670627411988165};
	 std::vector<double> band_energy = {e1_ays,e2_ays};
	

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
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=4;
	//g.read_ggmp("../example_graphs/ggm_sigma_no_tp/",ggm, max);
	g.read_ggmp("../example_graphs/ggm_sigma_nofock_notp/",ggm, max);
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	
	std::this_thread::sleep_for(std::chrono::milliseconds(1500)); 

	
	g.ggm_label(ggm,0);
	
bool hf = true; ///set to true if you are expecting a hatree or fock type interaction. True works for all type of graph so its set to default. Setting to false 
	//speeds up the process
	//construction 
mband mb(interaction,interaction_value,band_energy,hf);
int  min_ord = 2;		
int  max_ord = 3;

std::vector<double> Beta_ext_vec = {50};
std::vector<double> Mfreq_ext_vec;
for (int i =0;i<100;i++){Mfreq_ext_vec.push_back(i);}



std::vector<mband::sampler_collector> sigma_ToSum;
std::vector<AmiGraph::graph_t>sigma_FromGraph;
auto   startTime = std::chrono::high_resolution_clock::now();
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
            mband::sampler_collector sigma_collector;
            mb.sigma_sampler(ggm[i][j].graph_vec[k], sigma_collector);
		    sigma_ToSum.push_back(sigma_collector);
			sigma_FromGraph.push_back(ggm[i][j].graph_vec[k]);
			if(!sigma_collector.fermionic_edge_species.empty()){
				mband::output_collector output_collector;
				mb.calculate_sampled_sigma(ggm[i][j].graph_vec[k], sigma_collector,output_collector,Beta_ext_vec,Mfreq_ext_vec);
				std::string filename = "2bandH2_o"+std::to_string(i)+"_g" +std::to_string(j)+"_n" +std::to_string(k)+".txt";
				mb.write_output(filename,output_collector,Beta_ext_vec,Mfreq_ext_vec);
				
			}
			
        }
    }
}



auto   endTime = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
std::cout<<"Final Duration is " << duration;
/*
std::vector<double> Beta_ext_vec = {50};
std::vector<double> Mfreq_ext_vec;
for (int i =0;i<1;i++){Mfreq_ext_vec.push_back(i);}

for (int i = 0; i < sigma_FromGraph.size();i++){
	if(!sigma_ToSum[i].fermionic_edge_species.empty()){
		mband::output_collector output_collector;
		mb.calculate_sampled_sigma(sigma_FromGraph[i], sigma_ToSum[i],output_collector,Beta_ext_vec,Mfreq_ext_vec);
	}
}
*/	
/*
AmiGraph::graph_t  gself0 = ggm[5][0].graph_vec[0];    
std::string outputfile0 = "sto-10_0" + std::to_string(ord) + "_g" + std::to_string(group) + "_v" + std::to_string(0) + ".txt";


AmiGraph::graph_t  gself1 = ggm[2][0].graph_vec[1];    
std::string outputfile1 = "sto-10_0" + std::to_string(ord) + "_g" + std::to_string(group) + "_v" + std::to_string(1) + ".txt";

AmiGraph::graph_t  gself2 = ggm[2][0].graph_vec[2];    
std::string outputfile2 = "sto-10_0" + std::to_string(ord) + "_g" + std::to_string(group) + "_v" + std::to_string(2) + ".txt";
	
AmiGraph::graph_t  gself3 = ggm[2][0].graph_vec[3];    
std::string outputfile3 = "sto-10_0" + std::to_string(ord) + "_g" + std::to_string(group) + "_v" + std::to_string(3) + ".txt";

    

mband::sampler_collector Collector0;
mband::sampler_collector Collector1;
mband::sampler_collector Collector2;
mband::sampler_collector Collector3;	



mband::output_collector collector0;
mband::output_collector collector1;
mband::output_collector collector2;
mband::output_collector collector3;
*/


//std::vector<int> line = {2,2};
//mb.sigma_sampler(gself0, Collector0);

//mb.molecular_solver(gself0, collector0, Beta_ext_vec, Mfreq_ext_vec);




//mb.write_output(outputfile0, collector0, Beta_ext_vec, Mfreq_ext_vec);

print1d(band_energy);








//////////////////////////////////////////////////////////////////////////////////////////////////
////Does the same thing, but broken in step by step processs. Meant as a test bed and debugging!	
//////////////////////////////////////////////////////////////////////////////////////////////////

/*
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
	double prefactor =  g.get_prefactor(gself,ord);
	
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
	std::cout << "energy V are \n" <<std::endl;
	for (int i = 0; i <energy_vector.size();i++){
		for (int j = 0; j <energy_vector[i].size();j++)
		{
			std::cout<<energy_vector[i][j];
		}
		std::cout<<std::endl;
	}
	///evaluation using AMI _-_-_-_-_-


int n = 2*ord-1; // arbitrary number of elements in the array
AmiBase::g_struct gs[n];
std::vector<AmiBase::g_struct> gs_vec;
for(int i = 0; i < n; i++) {
     
    gs[i] = {Epsilon[i], Alpha[i]};
	gs_vec.push_back(gs[i]);
}

//AmiBase::g_struct gs[3] = {{Epsilon[0],Alpha[0]},{Epsilon[1],Alpha[1]},{Epsilon[2],Alpha[2]} };

AmiBase::g_prod_t R0 =gs_vec;
AmiBase::S_t S_array;
AmiBase::P_t P_array;
AmiBase::R_t R_array;

double E_REG=0; 
int N_INT=ord;
  
AmiBase::ami_parms test_amiparms(N_INT, E_REG);
//ami.drop_matsubara_poles = false;
ami.construct(test_amiparms , R0 , R_array , P_array , S_array ); 

AmiBase:: frequency_t frequency;


///
mband::output_collector collector;

std::vector<double> beta_ext_vec= {50};
std::vector<double> mfreq_ext_vec = {0};



std::complex<double> final_result = {0,0};
//Self energy for  different band indexes are evaluated for a set matsubara frequencies and beta
int count = 1;
for (int x = 0; x <mfreq_ext_vec.size();x++){ 
	for (int y = 0;y < beta_ext_vec.size();y++) {
		for (int i= 0; i<ord;i++){  frequency.push_back(std::complex<double>(0,0));}
		frequency.push_back(std::complex<double>(0,(2*mfreq_ext_vec [x]+ 1)*M_PI/beta_ext_vec[y]));
		    for (int i = 0; i<energy_vector.size();i++){
				
				
				AmiBase::ami_vars external (energy_vector[i],frequency,beta_ext_vec[y]);
				std::complex < double > calc_result = prefactor*ami.evaluate(test_amiparms,R_array ,P_array,S_array,external)*mb.Umatch(interaction,interaction_value,interaction_species [i]);
				
				
				collector.result_vec.push_back(calc_result);
				collector.beta_vec.push_back(beta_ext_vec[y] );
				collector.mfreq_vec.push_back((2*mfreq_ext_vec[x]+1)*M_PI/beta_ext_vec[y]);
				collector.Uindex_vec.push_back(mb.interaction_index(interaction_species[i]));
				collector.extline_vec.push_back(external_line[i]);			
				std::cout<<count <<". " << calc_result <<std::endl;
				//std::cout<<"result is " <<  calc_result <<std::endl;
				final_result = final_result+  calc_result;
				count++;
				}
			}
		}
		

std::cout <<"Pre-factor used is " << prefactor <<std::endl;
std::cout<<"Final Result with default U_matrix provided is " << final_result <<"\n";
std::cout<<"Number of internal fermionic line is:" <<fermionic_edge.size()<< std::endl;
std::cout<<"different number of possible species arrangements are: " <<fermionic_edge_species.size() <<std::endl;
std::cout<<"total number possible U_abcd interaction printed above: " <<interaction_species.size() <<"\n \n";


//print2d(external_line);
std::cout<<" Mfreq result is \n\n";
print1d(collector.mfreq_vec);
std::cout<<" beta result is \n\n";
print1d(collector.beta_vec);
std::cout<<"resut is \n\n";
print1d(collector.result_vec);
std::cout<<"Uindex is \n";
print2d(collector.Uindex_vec);
std::cout<<"external line species is \n";
print2d(collector.extline_vec);

std::ofstream outFile(outputfile);
for (int i = 0; i <collector.mfreq_vec.size();i++){
	outFile << collector.beta_vec[i] <<" " << collector.mfreq_vec[i] <<" " <<collector.result_vec[i].real() << " " <<collector.result_vec[i].imag() <<" " <<collector.extline_vec[i][0] <<" " << collector.extline_vec[i][1] <<" ";
	for (int j = 0; j< collector.Uindex_vec[i].size(); j++){
		outFile<< collector.Uindex_vec[i][j] <<" ";		
	}
	outFile << std::endl;
}
outFile.close();
*/

}


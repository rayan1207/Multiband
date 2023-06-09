
#include "mini_ami.hpp"

int main(int argc, char** argv)
{ 


    AmiBase ami;

   
     std::vector<std::vector<int>> interaction = readFile("data.txt");
	 std::vector<double> interaction_value = readFile1("data",5);
	 std::vector<double> band_energy = readFile1("data_h.txt",2);
	 
	 std::cout<<interaction.size();
    
	
	
	/*
    double e1 = -0.5825365736653964;
	double e2 = 0.6670627411988165;
	
	 double e1_ays  = 0.5*e1 - 0.5*e2;
	 double e2_ays  = 0.5*e2 -0.5*e1;

	std::cout <<"energies are" << e1_ays <<" " <<e2_ays;
	 std::vector<vector<int>> interaction = {{1,1,1,1},{2,2,1,1},{2,1,2,1},{1,2,2,1},{2,1,1,2},{1,2,1,2},{1,1,2,2},{2,2,2,2}};
     std::vector<vector<int>> hubbard = {{2,2,1,1},{1,1,2,2}};
     
	 std::vector<double> interaction_value={0.69907,0.664234,0.1815454,0.1815454,0.1815454,0.1815454,0.664234,0.674536};
	 std::vector<double> band_energy = {-0.5825365736653964,0.6670627411988165};
	 //std::vector<double> band_energy = {e1_ays,e2_ays};
  */
////////////////////////////James loading  stuff/////////////////////////////////
		int seed=3;	
		std::cout<< "Constructing AmiGraph object using seed: "<<seed <<" " <<std::endl;
		AmiGraph g(AmiBase::Sigma, seed);	

	std::string infile("ext_vars.dat");
	
		std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
		NewAmiCalc::external_variable_list extern_list;

	// read column data into a file - you can get external variables however you would like	
	g.ami.read_external(infile, extern_list);	
	



	AmiGraph::gg_matrix_t ggm; 
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=3;
	//g.read_ggmp("../example_graphs/ggm_sigma_nofock_notp/",ggm, max);
	g.read_ggmp("../example_graphs/ggm_all/",ggm, max);
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	
	std::this_thread::sleep_for(std::chrono::milliseconds(1500)); 

	
	g.ggm_label(ggm,0);
	

	
/////////////////////////////////////loading and labelling ends////////////////////




///////////////////////////////constructing mband ///////////////////////////////////////
	
bool hf = false; ///set to true if you are expecting a hatree or fock type interaction. True works for all type of graph so its set to default. Setting to false 
	//speeds up the process
	//construction 
mband mb(interaction,interaction_value,band_energy,hf);
int  min_ord = 2;		
int  max_ord = 2;

std::vector<double> Beta_ext_vec = {50};
std::vector<double> Mfreq_ext_vec;
for (int i =0;i<100;i++){Mfreq_ext_vec.push_back(i);}
std::vector<int> line={1,1};
std::vector<int> line1={2,2};
auto   startTime = std::chrono::high_resolution_clock::now();
std::vector<mband::sampler_collector> sigma_ToSum;
std::vector<AmiGraph::graph_t>sigma_FromGraph;

std::cout<<"External parameters read are"<<std::endl;
for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_freq_[0]<<std::endl;
	
		}	

///////////////////////////////just sampling//////////////////////////////////////////
/*
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
            mband::sampler_collector sigma_collector;
            mb.sigma_sampler(ggm[i][j].graph_vec[k], sigma_collector);
			sigma_FromGraph.push_back(ggm[i][j].graph_vec[k]);
			std::cout<<"the size is " <<sigma_collector.fermionic_edge_species.size();
		}
	}
}
*/

mband::sampler_collector sigma_collector;
mb.sigma_sampler(ggm[3][1].graph_vec[0], sigma_collector);

/*
int MC_num = 2000000;
int lattice_type = 1;
mband::sampler_collector sigma_collector;
mb.sigma_sampler(ggm[3][0].graph_vec[0], sigma_collector);

std::complex<double> result = mb.lcalc_sampled_sigma(ggm[3][0].graph_vec[0], sigma_collector.Epsilon, sigma_collector.Alpha, sigma_collector.fermionic_edge_species[0],
    extern_list[0],  MC_num, lattice_type);
	
	
std::cout <<" FINAL RESULT IS" << result/static_cast<double>(MC_num)<<std::endl;
*/



//////////////////////////////////sampling ends////////////////////////


///////////////////////////////using sampler,collecting and evaluating///////////////////////////////////
/*
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
            mband::sampler_collector sigma_collector;
            mb.sigma_sampler(ggm[i][j].graph_vec[k], sigma_collector);
		    sigma_ToSum.push_back(sigma_collector);
			sigma_FromGraph.push_back(ggm[i][j].graph_vec[k]);
			if(!sigma_collector.fermionic_edge_species.empty()){
				mband::output_collector output_collector;
				mb.calculate_sampled_sigma_ext(ggm[i][j].graph_vec[k], sigma_collector,output_collector,Beta_ext_vec,Mfreq_ext_vec,line);
				std::string filename = "10bandH2_o"+std::to_string(i)+"_g" +std::to_string(j)+"_n" +std::to_string(k)+"11.txt";
				mb.write_output(filename,output_collector,Beta_ext_vec,Mfreq_ext_vec);				
			}		
        }
    }
}

*/





///////////////////////////////////Molecular stuff on hold???//////////////////data  collected already/////////////////////////
/*
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
			mband::output_collector output_collector;
        	mb.molecular_solver_ext(ggm[i][j].graph_vec[k], output_collector, Beta_ext_vec, Mfreq_ext_vec,line);
            std::string filename = "2band_o"+std::to_string(i)+"_g" +std::to_string(j)+"_n" +std::to_string(k)+"11.txt";
		    mb.write_output(filename,output_collector,Beta_ext_vec,Mfreq_ext_vec);				
			}		
        }
    }

	
for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
			mband::output_collector output_collector;
        	mb.molecular_solver_ext(ggm[i][j].graph_vec[k], output_collector, Beta_ext_vec, Mfreq_ext_vec,line1);
            std::string filename = "2band_o"+std::to_string(i)+"_g" +std::to_string(j)+"_n" +std::to_string(k)+"22.txt";
		    mb.write_output(filename,output_collector,Beta_ext_vec,Mfreq_ext_vec);				
			}		
        }
    }
*/
//////////////////////////for a single graph///////////////
/*
mband::output_collector output_collector;
mb.molecular_solver_ext(ggm[3][0].graph_vec[10], output_collector, Beta_ext_vec, Mfreq_ext_vec,line);
std::string filename = "10band_o"+std::to_string(3)+"_g" +std::to_string(0)+"_n" +std::to_string(10)+"11.dat";
mb.write_output(filename,output_collector,Beta_ext_vec,Mfreq_ext_vec);

mband::output_collector output_collector1;
mb.molecular_solver_ext(ggm[3][0].graph_vec[11], output_collector1, Beta_ext_vec, Mfreq_ext_vec,line);
std::string filename1 = "10band_o"+std::to_string(3)+"_g" +std::to_string(0)+"_n" +std::to_string(11)+"11.dat";
mb.write_output(filename1,output_collector1,Beta_ext_vec,Mfreq_ext_vec);
*/
//////molecular stuff ends


auto   endTime = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
std::cout<<"Final Duration is " << duration;



}


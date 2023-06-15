
#include "mini_ami.hpp"

int main(int argc, char** argv)
{ 

AmiBase ami;
mband::params_param params;
params_loader("params.txt", params);
int seed =3;
AmiGraph g(AmiBase::Sigma, seed);
AmiGraph::gg_matrix_t ggm;	
NewAmiCalc::external_variable_list extern_list;
std::vector<std::vector<int>> interaction;
std::vector<double> interaction_value;
std::vector<double> band_energy;

MPI_Init(&argc, &argv);

int numProcesses, rank;
MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if (params.molecular==0&& params.lattice_type == 1){
    interaction = readFile("Hubbard_U.txt");
	interaction_value = readFile1("Hubbard_U.txt",5);
	std::vector<double> band_energy = {0,0};
    std::cout<< "Constructing AmiGraph object using seed: "<<seed <<" " <<std::endl;
	std::string infile("ext_vars.dat");
	std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
	g.ami.read_external(infile, extern_list);	
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=4;
	g.read_ggmp("../example_graphs/ggm_all/",ggm, max);
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	g.ggm_label(ggm,0); 
    std::cout<<"External parameters read are"<<std::endl;
for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_freq_[0]<<std::endl;
	}

	
}


else if (params.molecular==0&& params.lattice_type ==2){
     interaction = readFile("Hubbard_U.txt");
	 interaction_value = readFile1("Hubbard_U.txt",5);
	 std::vector<double> band_energy = {0,0};

    std::cout<< "Constructing AmiGraph object using seed: "<<seed <<" " <<std::endl;	
	std::string infile("ext_vars.dat");
	std::cout<<"Reading external parameters from ext_vars.dat"<<std::endl;
	g.ami.read_external(infile, extern_list);	
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	int max=4;
	g.read_ggmp("../example_graphs/ggm_sigma_no_tp/",ggm, max);
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	g.ggm_label(ggm,0);    
	
	std::cout<<"External parameters read are"<<std::endl;
    for(int i=0; i<extern_list.size();i++){
		std::cout<<extern_list[i].BETA_<<" "<<extern_list[i].MU_<<" "<< extern_list[i].H_<<" "<<extern_list[i].KDIM_<<" "<<extern_list[i].external_k_list_[0][0]<<" "<<extern_list[i].external_k_list_[0][1] <<" "<<extern_list[i].external_freq_[0]<<std::endl;
	
		}	

}
else {std::cout<<"Molecular stuff";}	



if (params.molecular == 0){
std::cout<<"beginning lattice calculations"<<std::endl;
bool hf = params.hatree_fock; ///set to true if you are expecting a hatree or fock type interaction. True works for all type of graph so its set to default. Setting to false 
	//speeds up the process
	//construction 
mband mb(interaction,interaction_value,band_energy,hf);	
int  min_ord = params.min_ord;		
int  max_ord = params.max_ord;
auto   startTime = std::chrono::high_resolution_clock::now();
std::vector<mband::sampler_collector> sigma_ToSum;
std::vector<AmiGraph::graph_t>sigma_FromGraph;

for (int i = min_ord; i < max_ord+1; ++i) {
    for (int j = 0; j < ggm[i].size(); ++j) {
        for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
            mband::sampler_collector sigma_collector;
            mb.sigma_sampler(ggm[i][j].graph_vec[k], sigma_collector);
			sigma_ToSum.push_back(sigma_collector);
			sigma_FromGraph.push_back(ggm[i][j].graph_vec[k]);
		}
	}
}

int MC_num = params.MC_num;
int lattice_type = params.lattice_type;

// Set the number of samples for each process
int samplesPerProcess = MC_num / numProcesses;
int remainingSamples = MC_num % numProcesses;

    // Adjust samplesPerProcess for the last process if numSamples is not divisible by numProcesses
if (rank == numProcesses - 1) {
        samplesPerProcess += remainingSamples;
    }
 std::complex<double> cmpzero(0.0, 0.0);

std::vector<std::vector<std::vector<std::complex<double>>>> localSums(sigma_ToSum.size(), std::vector<std::vector<std::complex<double>>>(extern_list.size()));
std::vector<std::vector<std::vector<std::complex<double>>>> localSumSquared(sigma_ToSum.size(), std::vector<std::vector<std::complex<double>>>(extern_list.size()));
std::vector<std::vector<std::vector<std::complex<double>>>> globalSums(sigma_ToSum.size(), std::vector<std::vector<std::complex<double>>>(extern_list.size()));
std::vector<std::vector<std::vector<std::complex<double>>>> globalSumSquared(sigma_ToSum.size(), std::vector<std::vector<std::complex<double>>>(extern_list.size()));



for (int i = 0; i < sigma_ToSum.size(); i++) {
    for (int j = 0; j < extern_list.size(); j++) {
        int speciesSize = sigma_ToSum[i].fermionic_edge_species.size();
        globalSums[i][j].resize(speciesSize, std::complex<double>(0.0, 0.0));
        localSums[i][j].resize(speciesSize, std::complex<double>(0.0, 0.0));
		globalSumSquared[i][j].resize(speciesSize, std::complex<double>(0.0, 0.0));
        localSumSquared[i][j].resize(speciesSize, std::complex<double>(0.0, 0.0));
        if (params.in == -1 && params.out== -1){
			for (int k = 0; k < speciesSize; k++) {		   
            localSums[i][j][k] = mb.lcalc_sampled_sigma(sigma_ToSum[i].graph, sigma_ToSum[i].Epsilon, sigma_ToSum[i].Alpha, sigma_ToSum[i].fermionic_edge_species[k],
                extern_list[j], samplesPerProcess, params).first;
			localSumSquared[i][j][k] = mb.lcalc_sampled_sigma(sigma_ToSum[i].graph, sigma_ToSum[i].Epsilon, sigma_ToSum[i].Alpha, sigma_ToSum[i].fermionic_edge_species[k],
                extern_list[j], samplesPerProcess, params).second;
			}		
		}
		else{
			for (int k = 0; k < speciesSize; k++) {
				if (sigma_ToSum[i].external_line[k][0] == params.in && sigma_ToSum[i].external_line[k][1]==params.out){
				localSums[i][j][k] = mb.lcalc_sampled_sigma(sigma_ToSum[i].graph, sigma_ToSum[i].Epsilon, sigma_ToSum[i].Alpha, sigma_ToSum[i].fermionic_edge_species[k],
					extern_list[j], samplesPerProcess, params).first;
				localSumSquared[i][j][k] = mb.lcalc_sampled_sigma(sigma_ToSum[i].graph, sigma_ToSum[i].Epsilon, sigma_ToSum[i].Alpha, sigma_ToSum[i].fermionic_edge_species[k],
					extern_list[j], samplesPerProcess, params).second;
				}
			}
		}
        
        MPI_Reduce(localSums[i][j].data(), globalSums[i][j].data(), speciesSize * 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(localSumSquared[i][j].data(), globalSumSquared[i][j].data(), speciesSize * 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
}

double U = 3;

if (rank == 0) {
std::cout << "molecular: " << params.molecular << std::endl;
std::cout << "E_reg: " << params.E_reg << std::endl;
std::cout << "lattice_type: " << params.lattice_type << std::endl;
std::cout << "set_precision: " << params.set_precision << std::endl;
std::cout << "hatree_fock: " << std::boolalpha << params.hatree_fock << std::endl;
std::cout << "tp: " << params.tp << std::endl;
std::cout << "tperp_p: " << params.tperp_p << std::endl;
std::cout << "tperp: " << params.tperp << std::endl;
std::cout << "Mc_num" << params.MC_num << std::endl;
std::cout << "max: " << params.max_ord << std::endl;
std::cout << "min_ord" << params.min_ord << std::endl;
std::cout << "line1 " << params.in<<" " <<params.out<< std::endl;
std::ofstream outputFile("result_test.txt"); // Open the result file for writing

    if (outputFile.is_open()) {
        std::cout << "Number of processes: " << numProcesses << std::endl;
        outputFile << "#Ord "<<"Beta " << "Mu " <<"H " <<"Kx " <<"Ky " <<"Rfq "<<"Imfq " <<"in " <<"out " <<"Re " << "Re_err " <<"Imag " << "Imag_err " <<"U1 " <<"U2 " <<"U3 "<<"U4 " <<std::endl;
        for (int i = 0; i < sigma_ToSum.size(); i++) {
            for (int j = 0; j < extern_list.size(); j++) {
                int speciesSize = sigma_ToSum[i].fermionic_edge_species.size();
                for (int k = 0; k < speciesSize; k++) {
					double Ufactor = std::pow(U,g.graph_order(sigma_ToSum[i].graph));
					
                    std::complex<double> integral =globalSums[i][j][k] / static_cast<double>(MC_num);
			        std::complex<double> integral_sq = globalSumSquared[i][j][k] ;
					double integral_err_real =  mb.Stdev(integral_sq.real(), integral.real(),MC_num );
					double integral_err_imag =  mb.Stdev(integral_sq.imag(), integral.imag(),MC_num );
					
					

                    std::cout << g.graph_order(sigma_ToSum[i].graph) << " ";
                    outputFile << g.graph_order(sigma_ToSum[i].graph) << " ";

                    std::cout << extern_list[j].BETA_ << " " << extern_list[j].MU_.real() << " " << extern_list[j].H_ << " "
                              << extern_list[j].external_k_list_[0][0] << " " << extern_list[j].external_k_list_[0][1]
                              << " " << extern_list[j].external_freq_[0].real() << " " << extern_list[j].external_freq_[0].imag()<<" ";
                    outputFile << extern_list[j].BETA_ << " " << extern_list[j].MU_.real() << " " << extern_list[j].H_ << " "
                               << extern_list[j].external_k_list_[0][0] << " " << extern_list[j].external_k_list_[0][1]
                               << " " << extern_list[j].external_freq_[0].real() << " " << extern_list[j].external_freq_[0].imag()<<" ";

                    std::cout << sigma_ToSum[i].external_line[k][0] << " " << sigma_ToSum[i].external_line[k][1] << " ";
                    outputFile << sigma_ToSum[i].external_line[k][0] << " " << sigma_ToSum[i].external_line[k][1] << " ";

                    std::cout << integral.real() << " "<<integral_err_real<<" "<<integral.imag() << " "<<integral_err_imag<<" " ;
                    outputFile <<integral.real() << " "<<integral_err_real<<" "<< integral.imag() << " "<<integral_err_imag<<" " ;

                    for (int l = 0; l < sigma_ToSum[i].Uindex[k].size(); l++) {
                        std::cout << sigma_ToSum[i].Uindex[k][l] << " ";
                        outputFile << sigma_ToSum[i].Uindex[k][l] << " ";
                    }
					int remainingElements = 4 - sigma_ToSum[i].Uindex[k].size();
						for (int l = 0; l < remainingElements; l++) {
							std::cout << -1 << " ";
							outputFile << -1 << " ";
														}

                    std::cout << std::endl;
                    outputFile << std::endl;
                }
            }
        }

        outputFile.close(); // Close the result file
	}
	
	else {
        std::cout << "Error: Failed to open the result file." << std::endl;
    }

	}
MPI_Finalize();
auto   endTime = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
std::cout<<"Final Duration is " << duration;
}
else{std::cout<<"molecular stuff";}


//////////////////////////////////sampling ends////////////////////////


///////////////////////////////using sampler,collecting and evaluating///////////////////////////////////
/*
std::vector<double> Beta_ext_vec = {50};
std::vector<double> Mfreq_ext_vec;
for (int i =0;i<100;i++){Mfreq_ext_vec.push_back(i);}
std::vector<int> line={1,1};
std::vector<int> line1={2,2};
*/
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






}


#include "mini_ami.hpp"

std::vector<std::complex < double >> result_collector;
std::vector<double> beta_collector;
std::vector<double> mfreq_collector;
std::vector<std::vector<int>> ext_line_collector;
std::vector<std::vector<int>>  Uindex_collector;
AmiGraph gg(AmiBase::Sigma, 0);
AmiBase Ami;


void mband::molecular_solver( AmiGraph::graph_t &gself, mband::output_collector& collector,std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec ){ 
	AmiGraph g(AmiBase::Sigma, 0);
	AmiBase ami;
	int ord = gg.graph_order(gself);
	AmiGraph::edge_vector_t fermionic_edge;
	std::vector<std::vector<int>> fermionic_edge_species;
	std::vector<std::vector<std::vector<int>>> interaction_species ;
	std::vector<std::vector<int>> external_line;
	std::vector<AmiBase::epsilon_t> Epsilon;
	std::vector<AmiBase::alpha_t> Alpha;
	mband::solve_multiband(gself,fermionic_edge,fermionic_edge_species,interaction_species ,external_line);
	mband::generate_eps_alpha(gself,fermionic_edge,Epsilon,Alpha);
	double prefactor =  g.get_prefactor(gself,ord);
	
	std::vector<std::vector<double>> speciesToEnergy = mband::band_to_hab(fermionic_edge_species);
	std::vector<std::vector<std::complex<double> >> energy_vector;
	for (auto vec: speciesToEnergy){
		energy_vector.push_back(mband::generate_ept(Epsilon, vec));

	}
	
	
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

	int n = 2*ord-1; //number of fermionic lines
	AmiBase::g_struct gs[n];
	std::vector<AmiBase::g_struct> gs_vec;
	for(int i = 0; i < n; i++) {
		 
		gs[i] = {Epsilon[i], Alpha[i]};
		gs_vec.push_back(gs[i]);
	}
	AmiBase::g_prod_t R0 =gs_vec;
	AmiBase::S_t S_array;
	AmiBase::P_t P_array;
	AmiBase::R_t R_array;

	double E_REG=0; 
	int N_INT=ord;
	  
	AmiBase::ami_parms test_amiparms(N_INT, E_REG);
	Ami.construct(test_amiparms , R0 , R_array , P_array , S_array ); 
	AmiBase::frequency_t frequency;
	
	
	std::complex<double> final_result = {0,0};

	int count = 1;
	for (int x = 0; x <mfreq_ext_vec.size();x++){ 
		for (int y = 0;y < beta_ext_vec.size();y++) {
			for (int i= 0; i<ord;i++){  frequency.push_back(std::complex<double>(0,0));}
			frequency.push_back(std::complex<double>(0,(2*mfreq_ext_vec [x]+ 1)*M_PI/beta_ext_vec[y]));
				for (int i = 0; i<energy_vector.size();i++){
					
					
					AmiBase::ami_vars external (energy_vector[i],frequency,beta_ext_vec[y]);
					std::complex < double > calc_result = prefactor*Ami.evaluate(test_amiparms,R_array ,P_array,S_array,external)*mband::Umatch(mband::interaction_legs,mband::int_values,interaction_species[i]);					
					collector.result_vec.push_back(calc_result);
					collector.beta_vec.push_back(beta_ext_vec[y] );
					collector.mfreq_vec.push_back((2*mfreq_ext_vec[x]+1)*M_PI/beta_ext_vec[y]);
					collector.Uindex_vec.push_back(mband::interaction_index(interaction_species[i]));
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
}

void mband::write_output(std::string outputfile,mband::output_collector& collector,std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec){ 
	std::ofstream outFile(outputfile);
	for (int i = 0; i <collector.mfreq_vec.size();i++){
		outFile << collector.beta_vec[i] <<" " << collector.mfreq_vec[i] <<" " <<collector.result_vec[i].real() << " " <<collector.result_vec[i].imag() <<" " <<collector.extline_vec[i][0] <<" " << collector.extline_vec[i][1] <<" ";
		for (int j = 0; j< collector.Uindex_vec[i].size(); j++){
			outFile<< collector.Uindex_vec[i][j] <<" ";		
		}
		outFile << std::endl;
	}
	outFile.close();
}
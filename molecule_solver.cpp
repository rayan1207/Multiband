#include "mini_ami.hpp"

std::vector<std::complex < double >> result_collector;
std::vector<double> beta_collector;
std::vector<double> mfreq_collector;
std::vector<std::vector<int>> ext_line_collector;
std::vector<std::vector<int>>  Uindex_collector;




void mband::molecular_solver( AmiGraph::graph_t &gself, mband::output_collector& collector,std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec,std::string outputfile ){ 
AmiGraph g(AmiBase::Sigma, 0);
	AmiBase ami;
	int ord = g.graph_order(gself);
	AmiGraph::edge_vector_t fermionic_edge;
	std::vector<std::vector<int>> fermionic_edge_species;
	std::vector<std::vector<std::vector<int>>> interaction_species ;
	std::vector<std::vector<int>> external_line;
	std::vector<AmiBase::epsilon_t> Epsilon;
	std::vector<AmiBase::alpha_t> Alpha;
	std::vector<std::vector<int>> bosonic_Alpha;
	mband::solve_multiband(gself,fermionic_edge,fermionic_edge_species,interaction_species ,bosonic_Alpha,external_line);
	mband::generate_eps_alpha(gself,fermionic_edge,Epsilon,Alpha);
	
	//double  loop_ord =std::pow( -1, (double) (g.count_fermi_loops(gself) + ord));
	
	double prefactor =  g.get_prefactor(gself,ord)*std::pow( 2, (double) g.count_fermi_loops(gself));
	
	std::vector<std::vector<double>> speciesToEnergy = mband::band_to_hab(fermionic_edge_species);
	
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
	double E_REG=0; 
	int N_INT=ord;  
	AmiBase::ami_parms test_amiparms(N_INT, E_REG);
	AmiBase::frequency_t frequency;
	std::complex<double> final_result = {0,0};
	std::vector<double> U;
	for (auto  i : interaction_species){
	 U.push_back(mband::Umatch(mband::interaction_legs,mband::int_values,i));
	}
	int count = 1;
    std::ofstream outFile(outputfile);
	for (int y = 0;y < beta_ext_vec.size();y++){ 
		for (int x = 0; x <mfreq_ext_vec.size();x++) {
			for (int i= 0; i<ord;i++){  frequency.push_back(std::complex<double>(0,0));}
			frequency.push_back(std::complex<double>(0,(2*mfreq_ext_vec [x]+ 1)*M_PI/beta_ext_vec[y]));
				for (int i = 0; i<speciesToEnergy.size();i++){					
					std::vector<std::complex<double>> energy_t = mband::generate_ept(Epsilon, speciesToEnergy[i]);
					AmiBase::ami_vars external (energy_t,frequency,beta_ext_vec[y]);
					AmiBase::g_struct gs[n];
						
					std::vector<AmiBase::g_struct> gs_vec;
					std::vector<AmiBase::epsilon_t> updated_Epsilon = mband::updateEpsilon(Epsilon,speciesToEnergy[i]);
	
					for(int i = 0; i < n; i++) {
		 
						gs[i] = {updated_Epsilon[i], Alpha[i]};
						gs_vec.push_back(gs[i]);
						}
					AmiBase::S_t S_array;
					AmiBase::P_t P_array;
					AmiBase::R_t R_array;
					AmiBase::g_prod_t R0 =gs_vec;
					ami.precision_cutoff=0;
					ami.drop_bosonic_diverge =true;
					ami.drop_matsubara_poles = false;
					ami.construct(test_amiparms , R0 , R_array , P_array , S_array ); 
					
					std::complex < double > calc_result = prefactor*ami.evaluate_otf(test_amiparms,R_array ,P_array,S_array,external)*U[i];					
					outFile << beta_ext_vec[y] <<" " <<(2*mfreq_ext_vec [x]+ 1)*M_PI/beta_ext_vec[y]<<" " <<calc_result.real() << " " <<calc_result.imag() <<" " <<external_line[i][0]<<" " << external_line[i][1] <<std::endl;
					 /*
					 if (std::isinf(calc_result.real()) || std::isnan(calc_result.real()) ||
						std::isinf(calc_result.imag()) || std::isnan(calc_result.imag())) {
						calc_result = std::complex<double>(0.0, 0.0);  // Set to zero
								}

					*/
					/*
					std::cout<< "result is " <<  calc_result <<std::endl;
					collector.result_vec.push_back(calc_result);
					collector.beta_vec.push_back(beta_ext_vec[y] );
					collector.mfreq_vec.push_back((2*mfreq_ext_vec[x]+1)*M_PI/beta_ext_vec[y]);
					//collector.Uindex_vec.push_back(mband::interaction_index(interaction_species[i]));
					std::vector<int> pp = {0,0};
					collector.Uindex_vec.push_back(pp);
					collector.extline_vec.push_back(external_line[i]);			

					final_result = final_result+  calc_result;
					count++;
					*/
						
					}
					frequency.clear();
				}
			}
			
    outFile.close();
	std::cout <<"Pre-factor used is " << prefactor <<std::endl;
	
	std::cout<<"Final Result with default U_matrix provided is " << final_result <<"\n";
	std::cout<<"Number of internal fermionic line is:" <<fermionic_edge.size()<< std::endl;
	std::cout<<"different number of possible species arrangements are: " <<fermionic_edge_species.size() <<std::endl;
	std::cout<<"total number possible U_abcd interaction printed above: " <<interaction_species.size() <<"\n \n";
	
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



void mband::molecular_solver_ext( AmiGraph::graph_t &gself, mband::output_collector& collector,std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec, std::vector<int> line,std::string outputfile ){ 
	AmiGraph g(AmiBase::Sigma, 0);
	AmiBase ami;
	int ord = g.graph_order(gself);
	AmiGraph::edge_vector_t fermionic_edge;
	std::vector<std::vector<int>> fermionic_edge_species;
	std::vector<std::vector<std::vector<int>>> interaction_species ;
	std::vector<std::vector<int>> external_line;
	std::vector<AmiBase::epsilon_t> Epsilon;
	std::vector<AmiBase::alpha_t> Alpha;
	std::vector<AmiBase::alpha_t> bosonic_Alpha;
	mband::solve_multiband(gself,fermionic_edge,fermionic_edge_species,interaction_species ,bosonic_Alpha,external_line);
	mband::generate_eps_alpha(gself,fermionic_edge,Epsilon,Alpha);
	
	//double  loop_ord =std::pow( -1, (double) (g.count_fermi_loops(gself) + ord));
	
	double prefactor =  g.get_prefactor(gself,ord)*std::pow( 2, (double) g.count_fermi_loops(gself));
	
	std::vector<std::vector<double>> speciesToEnergy = mband::band_to_hab(fermionic_edge_species);
	
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
	double E_REG=0; 
	int N_INT=ord;  
	AmiBase::ami_parms test_amiparms(N_INT, E_REG);
	AmiBase::frequency_t frequency;
	std::complex<double> final_result = {0,0};
	std::vector<double> U;
	for (auto  i : interaction_species){
	 U.push_back(mband::Umatch(mband::interaction_legs,mband::int_values,i));
	}
	int count = 1;
    std::ofstream outFile(outputfile);
	for (int y = 0;y < beta_ext_vec.size();y++){ 
		for (int x = 0; x <mfreq_ext_vec.size();x++) {
			for (int i= 0; i<ord;i++){  frequency.push_back(std::complex<double>(0,0));}
			frequency.push_back(std::complex<double>(0,(2*mfreq_ext_vec [x]+ 1)*M_PI/beta_ext_vec[y]));
				for (int i = 0; i<speciesToEnergy.size();i++){					
					if (external_line[i][0] ==line[0]& external_line[i][1] == line[1])
					{	
					std::vector<std::complex<double>> energy_t = mband::generate_ept(Epsilon, speciesToEnergy[i]);
					AmiBase::ami_vars external (energy_t,frequency,beta_ext_vec[y]);
					AmiBase::g_struct gs[n];
						
					std::vector<AmiBase::g_struct> gs_vec;
					std::vector<AmiBase::epsilon_t> updated_Epsilon = mband::updateEpsilon(Epsilon,speciesToEnergy[i]);
	
					for(int i = 0; i < n; i++) {
		 
						gs[i] = {updated_Epsilon[i], Alpha[i]};
						gs_vec.push_back(gs[i]);
						}
					AmiBase::S_t S_array;
					AmiBase::P_t P_array;
					AmiBase::R_t R_array;
					AmiBase::g_prod_t R0 =gs_vec;
					ami.precision_cutoff=0;
					ami.drop_bosonic_diverge =true;
					ami.drop_matsubara_poles = false;
					ami.construct(test_amiparms , R0 , R_array , P_array , S_array ); 
					
					std::complex < double > calc_result = prefactor*ami.evaluate_otf(test_amiparms,R_array ,P_array,S_array,external)*U[i];					
					outFile << beta_ext_vec[y] <<" " <<(2*mfreq_ext_vec [x]+ 1)*M_PI/beta_ext_vec[y]<<" " <<calc_result.real() << " " <<calc_result.imag() <<" " <<line[0] <<" " << line[1] <<std::endl;
					 /*
					 if (std::isinf(calc_result.real()) || std::isnan(calc_result.real()) ||
						std::isinf(calc_result.imag()) || std::isnan(calc_result.imag())) {
						calc_result = std::complex<double>(0.0, 0.0);  // Set to zero
								}

					*/
					/*
					std::cout<< "result is " <<  calc_result <<std::endl;
					collector.result_vec.push_back(calc_result);
					collector.beta_vec.push_back(beta_ext_vec[y] );
					collector.mfreq_vec.push_back((2*mfreq_ext_vec[x]+1)*M_PI/beta_ext_vec[y]);
					//collector.Uindex_vec.push_back(mband::interaction_index(interaction_species[i]));
					std::vector<int> pp = {0,0};
					collector.Uindex_vec.push_back(pp);
					collector.extline_vec.push_back(external_line[i]);			

					final_result = final_result+  calc_result;
					count++;
					*/
						}
					}
					frequency.clear();
				}
			}
			
    outFile.close();
	std::cout <<"Pre-factor used is " << prefactor <<std::endl;
	
	std::cout<<"Final Result with default U_matrix provided is " << final_result <<"\n";
	std::cout<<"Number of internal fermionic line is:" <<fermionic_edge.size()<< std::endl;
	std::cout<<"different number of possible species arrangements are: " <<fermionic_edge_species.size() <<std::endl;
	std::cout<<"total number possible U_abcd interaction printed above: " <<interaction_species.size() <<"\n \n";


}


void mband::sigma_sampler( AmiGraph::graph_t &gself, mband::sampler_collector& collector){ 
	AmiGraph g(AmiBase::Sigma, 0);
	AmiBase ami;
	int ord = g.graph_order(gself);

	auto   startTime = std::chrono::high_resolution_clock::now();
	mband::solve_multiband(gself,collector.fermionic_edge,collector.fermionic_edge_species,collector.interaction_species ,collector.bosonic_Alpha,collector.external_line);
	mband::generate_eps_alpha(gself,collector.fermionic_edge,collector.Epsilon,collector.Alpha);
	auto   endTime = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
	

	collector.graph = gself;
	
	for (auto interac: collector.interaction_species){
	collector.Uindex.push_back(interaction_index(interac));
	}
	/*
	for (int i=0;i <collector.fermionic_edge.size();i++){	
		g.print_edge_info(collector.fermionic_edge[i],gself);
		std::cout<<"epsilon we generate for AMI is ";
		print1d(collector.Epsilon[i]);
		std::cout<<" alpha is " ;
		print1d(collector.Alpha[i]);
		std::cout<<" all possible band indexes are: ";
		for (int j = 0; j < collector.fermionic_edge_species.size(); j++){
			std::cout<< collector.fermionic_edge_species[j][i] <<" ";
		}
	std::cout<<std::endl <<"\n";	
	}
	*/
   std::cout << " Number fermionic edge species" << std::endl;
   std::cout << "Execution time: " << duration << " milliseconds" << std::endl;
}

void mband::calculate_sampled_sigma(AmiGraph::graph_t &gself, mband::sampler_collector& samp_collector,  mband::output_collector& out_collector, std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec ){
   AmiGraph g(AmiBase::Sigma, 0);
   AmiBase ami;
   int ord = g.graph_order(gself);
   double prefactor =  g.get_prefactor(gself,ord)*std::pow( 2, (double) g.count_fermi_loops(gself));
   
   
   std::vector<std::vector<double>> speciesToEnergy = mband::band_to_hab(samp_collector.fermionic_edge_species);
   std::vector<std::vector<std::complex<double> >> energy_vector;
	for (auto vec: speciesToEnergy){
		energy_vector.push_back(mband::generate_ept(samp_collector.Epsilon, vec));
	}
	
	
	int n = 2*ord-1; //number of fermionic lines
	AmiBase::g_struct gs[n];
	std::vector<AmiBase::g_struct> gs_vec;
	for(int i = 0; i < n; i++) {
		 
		gs[i] = {samp_collector.Epsilon[i], samp_collector.Alpha[i]};
		gs_vec.push_back(gs[i]);
	}
	AmiBase::g_prod_t R0 =gs_vec;
	AmiBase::S_t S_array;
	AmiBase::P_t P_array;
	AmiBase::R_t R_array;

	double E_REG=1e-8; 
	int N_INT=ord;
	  
	AmiBase::ami_parms test_amiparms(N_INT, E_REG);
	ami.construct(test_amiparms , R0 , R_array , P_array , S_array ); 
	AmiBase::frequency_t frequency;
	
	
	std::complex<double> final_result = {0,0};


	
	for (int y = 0;y < beta_ext_vec.size();y++){ 
		for (int x = 0; x <mfreq_ext_vec.size();x++) {
			for (int i= 0; i<ord;i++){  frequency.push_back(std::complex<double>(0,0));}
			frequency.push_back(std::complex<double>(0,(2*mfreq_ext_vec [x]+ 1)*M_PI/beta_ext_vec[y]));
				for (int i = 0; i<energy_vector.size();i++){
					AmiBase::ami_vars external (energy_vector[i],frequency,beta_ext_vec[y]);
					std::complex < double > calc_result = prefactor*ami.evaluate(test_amiparms,R_array ,P_array,S_array,external)*mband::Umatch(mband::interaction_legs,
					mband::int_values,samp_collector.interaction_species[i]);					
					out_collector.result_vec.push_back(calc_result);
					out_collector.beta_vec.push_back(beta_ext_vec[y] );
					out_collector.mfreq_vec.push_back((2*mfreq_ext_vec[x]+1)*M_PI/beta_ext_vec[y]);
					out_collector.Uindex_vec.push_back(mband::interaction_index(samp_collector.interaction_species[i]));
					out_collector.extline_vec.push_back(samp_collector.external_line[i]);			
					std::cout<<"result is " <<  calc_result <<std::endl;
					final_result = final_result+  calc_result;
					}
					frequency.clear();
				}
			}

	std::cout <<"Pre-factor used is " << prefactor <<std::endl;
	std::cout<<"Final Result with default U_matrix provided is " << final_result <<"\n";
	std::cout<<"Number of internal fermionic line is:" <<samp_collector.fermionic_edge.size()<< std::endl;
	std::cout<<"different number of possible species arrangements are: " <<samp_collector.fermionic_edge_species.size() <<std::endl;
	std::cout<<"total number possible U_abcd interaction printed above: " <<samp_collector.interaction_species.size() <<"\n \n";
	

}	


void mband::calculate_sampled_sigma_ext(AmiGraph::graph_t &gself, mband::sampler_collector& samp_collector,  mband::output_collector& out_collector, std::vector<double> beta_ext_vec,std::vector<double> mfreq_ext_vec,std::vector<int> line ){
   AmiGraph g(AmiBase::Sigma, 0);
   AmiBase ami;
   int ord = g.graph_order(gself);
   double prefactor =  g.get_prefactor(gself,ord)*std::pow( 2, (double) g.count_fermi_loops(gself));
   
   
   std::vector<std::vector<double>> speciesToEnergy = mband::band_to_hab(samp_collector.fermionic_edge_species);
   std::vector<std::vector<std::complex<double> >> energy_vector;
	for (auto vec: speciesToEnergy){
		energy_vector.push_back(mband::generate_ept(samp_collector.Epsilon, vec));
	}
	
	
	int n = 2*ord-1; //number of fermionic lines
	AmiBase::g_struct gs[n];
	std::vector<AmiBase::g_struct> gs_vec;
	for(int i = 0; i < n; i++) {
		 
		gs[i] = {samp_collector.Epsilon[i], samp_collector.Alpha[i]};
		gs_vec.push_back(gs[i]);
	}
	AmiBase::g_prod_t R0 =gs_vec;
	AmiBase::S_t S_array;
	AmiBase::P_t P_array;
	AmiBase::R_t R_array;

	double E_REG=1e-8; 
	int N_INT=ord;
	  
	AmiBase::ami_parms test_amiparms(N_INT, E_REG);
	ami.construct(test_amiparms , R0 , R_array , P_array , S_array ); 
	AmiBase::frequency_t frequency;
	ami.precision_cutoff =1e-8;
	
	std::complex<double> final_result = {0,0};


	
	for (int y = 0;y < beta_ext_vec.size();y++){ 
		for (int x = 0; x <mfreq_ext_vec.size();x++) {
			for (int i= 0; i<ord;i++){  frequency.push_back(std::complex<double>(0,0));}
			frequency.push_back(std::complex<double>(0,(2*mfreq_ext_vec [x]+ 1)*M_PI/beta_ext_vec[y]));
				for (int i = 0; i<energy_vector.size();i++){
					if (samp_collector.external_line[i][0] ==line[0] & samp_collector.external_line[i][1] == line[1]){
						AmiBase::ami_vars external (energy_vector[i],frequency,beta_ext_vec[y]);	
						std::complex < double > calc_result = prefactor*ami.evaluate(test_amiparms,R_array ,P_array,S_array,external)*mband::Umatch(mband::interaction_legs,
						mband::int_values,samp_collector.interaction_species[i]);					
						out_collector.result_vec.push_back(calc_result);
						out_collector.beta_vec.push_back(beta_ext_vec[y] );
						out_collector.mfreq_vec.push_back((2*mfreq_ext_vec[x]+1)*M_PI/beta_ext_vec[y]);
						std::vector<int> pp ={0,0};
						//out_collector.Uindex_vec.push_back(mband::interaction_index(samp_collector.interaction_species[i]));
						out_collector.Uindex_vec.push_back(pp);
						out_collector.extline_vec.push_back(samp_collector.external_line[i]);			
						std::cout<<"result is " <<  calc_result <<std::endl;
					final_result = final_result+  calc_result;}
					
					}
					frequency.clear();
				}
			}

	std::cout <<"Pre-factor used is " << prefactor <<std::endl;
	std::cout<<"Final Result with default U_matrix provided is " << final_result <<"\n";
	std::cout<<"Number of internal fermionic line is:" <<samp_collector.fermionic_edge.size()<< std::endl;
	std::cout<<"different number of possible species arrangements are: " <<samp_collector.fermionic_edge_species.size() <<std::endl;
	std::cout<<"total number possible U_abcd interaction printed above: " <<samp_collector.interaction_species.size() <<"\n \n";
}


/*
std::complex<double> mband::lcalc_sampled_sigma(AmiGraph::graph_t &gself, std::vector<AmiBase::epsilon_t>& Epsilon, std::vector<AmiBase::alpha_t>& Alpha, std::vector<int>& Species,
    NewAmiCalc::ext_vars& ext_params, int MC_num, int lattice_type) {
    AmiGraph g(AmiBase::Sigma, 0);
    AmiBase ami;
    int ord = g.graph_order(gself);
    double prefactor = g.get_prefactor(gself,ord);
    int n = 2 * ord - 1; // number of fermionic lines
    AmiBase::g_struct gs[n];
    std::vector<AmiBase::g_struct> gs_vec;
    for (int i = 0; i < n; i++) {
        gs[i] = {Epsilon[i], Alpha[i]};
        gs_vec.push_back(gs[i]);
    }

    AmiBase::g_prod_t R0 = gs_vec;
    AmiBase::S_t S_array;
    AmiBase::P_t P_array;
    AmiBase::R_t R_array;

    double E_REG = 0;
    int N_INT = ord;

    AmiBase::ami_parms test_amiparms(N_INT, E_REG);
    ami.construct(test_amiparms, R0, R_array, P_array, S_array);
    AmiBase::frequency_t frequency;
    for (int i = 0; i < ord; i++) {
        frequency.push_back(std::complex<double>(0, 0));
    }
    frequency.push_back(ext_params.external_freq_[0]);
     
    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    int kspace = Alpha[0].size() - 1;
    std::complex<double> localSum(0.0, 0.0);
    // Generate random samples and calculate the sum
    for (int i = 0; i < MC_num; i++) {
        std::vector<std::vector<double>> momenta;
        momenta.reserve(kspace);
        for (int j = 0;  j< kspace; j++) {
            double momentum1 = 2 * M_PI * distribution(engine);
            double momentum2 = 2 * M_PI * distribution(engine);
            momenta.push_back({momentum1, momentum2});
        }
        momenta.push_back(ext_params.external_k_list_[0]);
		
		
        std::vector<std::vector<double>> summed_momenta;
        summed_momenta.reserve(Alpha.size());
        for (const auto& alpha : Alpha) {
            double qx = 0;
            double qy = 0;
            for (int j = 0; j < alpha.size(); j++) {
                qx += static_cast<double>(alpha[j]) * momenta[j][0];
                qy += static_cast<double>(alpha[j]) * momenta[j][1];
            }
            summed_momenta.push_back({qx, qy});
        }
	
        std::vector<double> energy;
        energy.reserve(summed_momenta.size());
        for (int i = 0; i < summed_momenta.size(); i++) {
            if (lattice_type == 1) {  
                energy.push_back(mband::Hubbard_Energy(ext_params, summed_momenta[i], Species[i]));
            }
        }

        std::vector<std::complex<double>> energy_t = mband::generate_ept(Epsilon, energy);
	
        AmiBase::ami_vars external(energy_t, frequency, ext_params.BETA_);
        localSum += prefactor* ami.evaluate(test_amiparms, R_array, P_array, S_array, external);
    }


	return localSum;
}

*/


std::pair<std::complex<double>, std::complex<double>> mband::lcalc_sampled_sigma(AmiGraph::graph_t &gself, std::vector<AmiBase::epsilon_t>& Epsilon, std::vector<AmiBase::alpha_t>& Alpha,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<int> &Utype,
 std::vector<int>& Species,NewAmiCalc::ext_vars& ext_params,int MC_num,params_param& param) {
    AmiGraph g(AmiBase::Sigma, 0);
    AmiBase ami;
    int ord = g.graph_order(gself);
    double prefactor = g.get_prefactor(gself, ord);
    int n = 2 * ord - 1; // number of fermionic lines
    AmiBase::g_struct gs[n];
    std::vector<AmiBase::g_struct> gs_vec;
    for (int i = 0; i < n; i++) {
        gs[i] = {Epsilon[i], Alpha[i]};
        gs_vec.push_back(gs[i]);
    }

    AmiBase::g_prod_t R0 = gs_vec;
    AmiBase::S_t S_array;
    AmiBase::P_t P_array;
    AmiBase::R_t R_array;

    double E_REG = param.E_reg;
    int N_INT = ord;
    ami.precision_cutoff=param.set_precision;
    AmiBase::ami_parms test_amiparms(N_INT, E_REG);
    ami.construct(test_amiparms, R0, R_array, P_array, S_array);
    AmiBase::frequency_t frequency;
    for (int i = 0; i < ord; i++) {
        frequency.push_back(std::complex<double>(0, 0));
    }
    frequency.push_back(ext_params.external_freq_[0]);

    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    int kspace = Alpha[0].size() - 1;
    std::complex<double> localSum(0.0, 0.0);
    std::complex<double> localSumOfSquares(0.0, 0.0);
    // Generate random samples and calculate the sum
    for (int i = 0; i < MC_num; i++) {
        std::vector<std::vector<double>> momenta;
        momenta.reserve(kspace);
        for (int j = 0; j < kspace; j++) {
            double momentum1 = 2 * M_PI * distribution(engine);
            double momentum2 = 2 * M_PI * distribution(engine);
            momenta.push_back({momentum1, momentum2});
        }
        momenta.push_back(ext_params.external_k_list_[0]);

        std::vector<std::vector<double>> summed_momenta;
        summed_momenta.reserve(Alpha.size());
        for (const auto& alpha : Alpha) {
            double qx = 0;
            double qy = 0;
            for (int j = 0; j < alpha.size(); j++) {
                qx += static_cast<double>(alpha[j]) * momenta[j][0];
                qy += static_cast<double>(alpha[j]) * momenta[j][1];
            }
            summed_momenta.push_back({qx, qy});
        }
	double form_factor = 1;
    
	/////////////////////////form factor for exteded hubbard////////////////////
    /*
	if (param.lattice_type ==3){		
			std::vector<std::vector<int>> nonlocal_alpha= mband::find_non_local_bosonic_alpha(bosonic_Alpha, Utype);
			std::vector<std::vector<double>> V_momenta;
			V_momenta.reserve(nonlocal_alpha.size());
        if (!nonlocal_alpha.empty()){
			for (const auto& b_alpha : nonlocal_alpha) {
				double qx = 0;
				double qy = 0;
				for (int j = 0; j < b_alpha.size(); j++) {
					qx += static_cast<double>(b_alpha[j]) * momenta[j][0];
					qy += static_cast<double>(b_alpha[j]) * momenta[j][1];
				}
			   V_momenta.push_back({ qx, qy });
			}

		    
			for (auto vq : V_momenta) {
				form_factor = form_factor*mband::non_local_U_formfactor(vq);
				}	
		}		
	}
	*/
	if (param.lattice_type ==3){		
		std::vector<std::vector<double>> V_momenta;
		V_momenta.reserve(bosonic_Alpha.size());

		for (const auto& b_alpha : bosonic_Alpha) {
			double qx = 0;
			double qy = 0;
			for (int j = 0; j < b_alpha.size(); j++) {
				qx += static_cast<double>(b_alpha[j]) * momenta[j][0];
				qy += static_cast<double>(b_alpha[j]) * momenta[j][1];
			}
		   V_momenta.push_back({ qx, qy });
		}

		for (int i = 0; i<Utype.size();i++){
			if (Utype[i]==0 || Utype[i]==1){
				form_factor= form_factor*2*param.V*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1]));
			}
			else {
				form_factor= form_factor*(1.0+  2*param.V*(std::cos(V_momenta[i][0])+std::cos(V_momenta[i][1])));	
	
			}	
		}
	}
		

	/////////////////////////////////energy/////////////////////////////	
        std::vector<double> energy;
        energy.reserve(summed_momenta.size());
        for (int i = 0; i < summed_momenta.size(); i++) {
            if (param.lattice_type == 1 || param.lattice_type == 3 ) {
                energy.push_back(mband::Hubbard_Energy(ext_params, summed_momenta[i], Species[i],param));
            }
			if (param.lattice_type == 2 ) {
                energy.push_back(mband::Bilayer_Hubbard_Energy(ext_params, summed_momenta[i], Species[i],param));
            }
        }

        std::vector<std::complex<double>> energy_t = mband::generate_ept(Epsilon, energy);

        AmiBase::ami_vars external(energy_t, frequency, ext_params.BETA_);
        std::complex<double> result =form_factor* prefactor * ami.evaluate(test_amiparms, R_array, P_array, S_array, external);
        localSum += result;
        localSumOfSquares += std::complex<double> (std::pow(result.real(),2),std::pow(result.imag(),2)) ;
    }

    return std::make_pair(localSum, localSumOfSquares);
}
	

     


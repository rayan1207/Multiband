#include "mini_ami.hpp"


AmiGraph g(AmiBase::Sigma, 0);
AmiBase ami;

///constructor
mband::mband(std::vector<std::vector<int>> _interaction_legs, std::vector<double> _int_values,std::vector<double> _energy, bool _Hartee_fock)
:  interaction_legs(_interaction_legs), int_values(_int_values),energy(_energy),  Hartee_fock (_Hartee_fock) {};

///from the interaction lines, 
std::vector<int> mband::Hartee_fock_filter(AmiGraph::edge_vector_t &fermionic_edge){
	std::vector<int> equalIndices;
	if (fermionic_edge[0] == fermionic_edge[1]){
		equalIndices= {0,1};
		return equalIndices;
	}
	
	else if (fermionic_edge[2] == fermionic_edge[3]){
		equalIndices= {2,3};
		return equalIndices;
		
	}
	else if (fermionic_edge[0] == fermionic_edge[2]){
		equalIndices= {0,2};
		return equalIndices;
	}
	
	else if (fermionic_edge[0] == fermionic_edge[3]){
		equalIndices= {0,3};
		return equalIndices;
	}
	else if (fermionic_edge[1] == fermionic_edge[2]){
		equalIndices= {1,2};
		return equalIndices;
	}
	
	else if (fermionic_edge[1] == fermionic_edge[3]){
		equalIndices= {1,3};
		return equalIndices;
	}
	else {
		return  equalIndices;
	}
	
	
}


void mband::find_interaction(AmiGraph::graph_t &graph, AmiGraph::edge_vector_t &b_vector, std::vector<AmiGraph::edge_vector_t> &f_vector){
	//g.find_bosonic_edges(graph, b_vector);
	boost::graph_traits<AmiGraph::graph_t>::in_edge_iterator iei, iedge_end;
    boost::graph_traits<AmiGraph::graph_t>::out_edge_iterator oei, oedge_end;
	AmiGraph::edge_vector_t v;
	for (int i=0;i <b_vector.size();i++){
		for (boost::tie(iei, iedge_end) = in_edges(source(b_vector[i],graph), graph); iei != iedge_end; ++iei)
    {    	if( graph[*iei].g_struct_.stat_==AmiBase::Fermi){
	v.push_back(*iei);}
	}
          
    
	for (boost::tie(oei, oedge_end) = out_edges(source(b_vector[i],graph), graph); oei != oedge_end; ++oei)
    {   
        if( graph[*oei].g_struct_.stat_==AmiBase::Fermi){
		v.push_back(*oei);}
    }
	for (boost::tie(iei, iedge_end) = in_edges(target(b_vector[i],graph), graph); iei != iedge_end; ++iei)
    {    if( graph[*iei].g_struct_.stat_==AmiBase::Fermi){
	       v.push_back(*iei);}
    }
	
	for (boost::tie(oei, oedge_end) = out_edges(target(b_vector[i],graph), graph); oei != oedge_end; ++oei)
    {   
        if( graph[*oei].g_struct_.stat_==AmiBase::Fermi){
		   v.push_back(*oei);}
    }

	
	f_vector.push_back(v);
	
	v.clear();
	  
	}
	
}



void mband::assign_label(AmiGraph::graph_t &g1, AmiGraph::edge_vector_t edge,vector<int> vector){
///assign species value given a vector edge {1,2,3,4}
for (int i =0; i< vector.size(); i++){
	
	
	g1[edge[i]].g_struct_.species_ = vector[i];

}
/*
for (int n : vector){std::cout << n ;}
std::cout <<std::endl;
*/
}

vector<int> mband::generate_edge_species(AmiGraph::graph_t &g, AmiGraph::edge_vector_t edge){
	std::vector<int> vector;
	for (int i =0; i< edge.size(); i++){
	
	
	vector.push_back(g[edge[i]].g_struct_.species_) ;

}


 return vector;
	
}

void mband::reset_species(AmiGraph::graph_t &graph,std::vector<AmiGraph::edge_vector_t> int_vector){
	for (int i = 0; i <int_vector.size();i++){
		mband::assign_label(graph,int_vector[i],{0,0,0,0});
	}
}

std::vector<std::vector<int>> mband::findmatch(std::vector<int> v1,AmiGraph::edge_vector_t &inter_vec) {
	std::vector<std::vector<int>> matchingVectors;
	if (mband::Hartee_fock){
		 matchingVectors = mband::interaction_legs;
		mband::filter(matchingVectors,mband::Hartee_fock_filter(inter_vec));
	}
	else {
		
		 matchingVectors = mband::interaction_legs;
		
	}

    for (int j = 0; j < v1.size(); j++) {
        if (v1[j] != 0) {
            int matchCount = 0;
            for (int i = 0; i < matchingVectors.size(); i++) {
                if (matchingVectors[i][j] == v1[j]) {
                    matchingVectors[matchCount] = matchingVectors[i];
                    matchCount++;
                }
            }
            matchingVectors.resize(matchCount);
        }
    }
 return matchingVectors;
}


std::vector<int> mband::external_species(AmiGraph::graph_t &graph){
	AmiGraph::vertex_vector_t v;
	AmiGraph::edge_vector_t edges;
	g.find_external_vertices(graph,v,edges);
	//std::cout<<"{" << graph[edges[0]].g_struct_.species_ <<"," << graph[edges[1]].g_struct_.species_<<"}";
	return {graph[edges[0]].g_struct_.species_,graph[edges[1]].g_struct_.species_};
}

void mband::solve_multiband_4(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs){
    AmiGraph::edge_vector_t bvector;
	std::vector<AmiGraph::edge_vector_t> int_vector;
	g.find_bosonic_edges(graph,bvector);
	g.find_internal_fermionic_edges(graph,fermionic_edge);
	for (int i=0;i <bvector.size();i++){
		//g.print_edge_info(bvector[i],gself2);
		std::cout << "Edge = ("  << graph[source(bvector[i],graph)].index_<<"," << graph[target(bvector[i],graph)].index_ <<")" 
		<<std::endl;
	}
	 for (auto b: bvector){
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
	}
	mband::find_interaction(graph,bvector,int_vector);
	mband::print_interactions(graph,bvector, int_vector);
	std::vector<std::vector<int>> possible_species_1;
	std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);
	if (mband::Hartee_fock){
		possible_species_1 = findmatch(initial_species_1,int_vector[0]);
	}
	else{
		possible_species_1 = mband::interaction_legs;
	}
	
	std::cout<< " Possible external species o4 in curly bracket and internal species in () is shown below \n";
	if (!possible_species_1.empty()){
		for (int i = 0; i<possible_species_1.size();i++){
			
			std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);	
			mband::assign_label(graph,int_vector[0],possible_species_1[i]);
			std::vector<int> initial_species_2 = mband::generate_edge_species(graph, int_vector[1]);
			std::vector<std::vector<int>> possible_species_2 = mband::findmatch(initial_species_2,int_vector[1]);
		

			if (!possible_species_2.empty()){
				for (int j = 0; j < possible_species_2.size();j++){
				//mband::print_match(possible_species_2[j],2);
				mband::assign_label(graph,int_vector[1],initial_species_2);				
				mband::assign_label(graph,int_vector[1],possible_species_2[j]);
				std::vector<int> initial_species_3 = mband::generate_edge_species(graph, int_vector[2]);
				std::vector<std::vector<int>> possible_species_3 = mband::findmatch(initial_species_3,int_vector[2]);
				
					if (!possible_species_3.empty()){
						for (int k = 0; k < possible_species_3.size();k++){
							//mband::print_match(possible_species_3[k],3);
							mband::assign_label(graph,int_vector[2],initial_species_3);
							mband::assign_label(graph,int_vector[2],possible_species_3[k]);											
						//interaction_species.push_back({mband::interaction_legs[i],possible_species_2[j],possible_species_3[k]});
											
							std::vector<int> initial_species_4 = mband::generate_edge_species(graph, int_vector[3]);
							std::vector<std::vector<int>> possible_species_4 = mband::findmatch(initial_species_4,int_vector[3]);					
							if (!possible_species_4.empty()){
								for (int l = 0; l < possible_species_4.size();l++){
									//num++
								//mband::print_match(possible_species_4[l],4);
								mband::assign_label(graph,int_vector[3],initial_species_4);							
								mband::assign_label(graph,int_vector[3],possible_species_4[l]);
								interaction_species.push_back({possible_species_1[i],possible_species_2[j],possible_species_3[k],possible_species_4[l]});
								ext_legs.push_back(mband::external_species(graph));
								//std::cout<< fermionic_edge.size();
								//std::cout<<j << " printing internal fermionic band index" << std::endl;
								//std::cout<<"(";		
								std::vector<int> v;
								for (int x =0; x < fermionic_edge.size(); x++){
									v.push_back(graph[fermionic_edge[x]].g_struct_.species_);
									//std::cout<<graph[fermionic_edge[x]].g_struct_.species_;
					
				}
				fermionic_species.push_back(v);
				v.clear();
								//std::cout<<")"<<std::endl;
								}mband::assign_label(graph,int_vector[3],initial_species_4);
								
							}
							else { mband::assign_label(graph,int_vector[3],initial_species_4);}				
						}mband::assign_label(graph,int_vector[2],initial_species_3);
										
					}
					else{mband::assign_label(graph,int_vector[2],initial_species_3);}
				}
				mband::assign_label(graph,int_vector[1],initial_species_2);			
			}
			else {mband::assign_label(graph,int_vector[1],initial_species_2);}
			
		mband::reset_species(graph,int_vector);	
		}
	}
	else{std::cout<< "NO match found \n";}
mband::print_assigned_species(interaction_species);
}
	

	
void mband::solve_multiband_3(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs){
    AmiGraph::edge_vector_t bvector;
	std::vector<AmiGraph::edge_vector_t> int_vector;
	g.find_bosonic_edges(graph,bvector);
	g.find_internal_fermionic_edges(graph,fermionic_edge);
    for (auto b: bvector){
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
	}	
	for (int i=0;i <bvector.size();i++){
		//g.print_edge_info(bvector[i],gself2);
		std::cout << "Edge = ("  << graph[source(bvector[i],graph)].index_<<"," << graph[target(bvector[i],graph)].index_ <<")" 
		<<std::endl;
	}
	mband::find_interaction(graph,bvector,int_vector);
	mband::print_interactions(graph,bvector, int_vector);
	std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);
	std::vector<std::vector<int>> possible_species_1;
	if (mband::Hartee_fock){
		possible_species_1 = findmatch(initial_species_1,int_vector[0]);
	}
	else{
		possible_species_1 = mband::interaction_legs;
	}
	std::cout<< " Possible external species in curly bracket and internal species in () is shown below \n";
	if (!possible_species_1.empty()) {
		for (int i = 0; i<possible_species_1.size();i++){
			std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);	
			mband::assign_label(graph,int_vector[0],possible_species_1[i]);

			std::vector<int> initial_species_2 = mband::generate_edge_species(graph, int_vector[1]);

			std::vector<std::vector<int>> possible_species_2 = mband::findmatch(initial_species_2,int_vector[1]);
			if (!possible_species_2.empty()){
				for (int j = 0; j < possible_species_2.size();j++){
				mband::assign_label(graph,int_vector[1],initial_species_2);				
				mband::assign_label(graph,int_vector[1],possible_species_2[j]);
				std::vector<int> initial_species_3 = mband::generate_edge_species(graph, int_vector[2]);
				std::vector<std::vector<int>> possible_species_3 = mband::findmatch(initial_species_3, int_vector[2]);
					if (!possible_species_3.empty()){
						for (int k = 0; k < possible_species_3.size();k++){
							mband::assign_label(graph,int_vector[2],initial_species_3);
							mband::assign_label(graph,int_vector[2],possible_species_3[k]);											
							interaction_species.push_back({possible_species_1[i],possible_species_2[j],possible_species_3[k]});
							ext_legs.push_back(mband::external_species(graph));
							//std::cout<< fermionic_edge.size();
							//std::cout<<j << " printing internal fermionic band index" << std::endl;
							//std::cout<<"(";
							std::vector<int> v;
				for (int x =0; x < fermionic_edge.size(); x++){
					v.push_back(graph[fermionic_edge[x]].g_struct_.species_);
					//std::cout<<graph[fermionic_edge[x]].g_struct_.species_ ;
				}
				
				fermionic_species.push_back(v);
				v.clear();
				//std::cout<<")"<<std::endl;
															
													
						}mband::assign_label(graph,int_vector[2],initial_species_3);
										
					}
					else{mband::assign_label(graph,int_vector[2],initial_species_3);}
				}
				mband::assign_label(graph,int_vector[1],initial_species_2);
				
			}
			else {mband::assign_label(graph,int_vector[1],initial_species_2);}
			
		mband::reset_species(graph,int_vector);
			}
    }
	else {std::cout << " No match found \n";}
mband::print_assigned_species(interaction_species);  
}







void mband::solve_multiband_2(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs){
    AmiGraph::edge_vector_t bvector;
	std::vector<AmiGraph::edge_vector_t> int_vector;
	g.find_bosonic_edges(graph,bvector);
	g.find_internal_fermionic_edges(graph,fermionic_edge);
	  for (auto b: bvector){
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
	}	
	for (int i=0;i <bvector.size();i++){
		//g.print_edge_info(bvector[i],gself2);
		std::cout << "Edge = ("  << graph[source(bvector[i],graph)].index_<<"," << graph[target(bvector[i],graph)].index_ <<")" 
		<<std::endl;
	}
	for (auto b: bvector){
		print1d(graph[b].g_struct_.alpha_);
	}
	mband::find_interaction(graph,bvector,int_vector);
	mband::print_interactions(graph,bvector, int_vector);
	std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);
	std::vector<std::vector<int>> possible_species_1;
	if (mband::Hartee_fock){
		possible_species_1 = findmatch(initial_species_1,int_vector[0]);
	}
	else{
		possible_species_1 = mband::interaction_legs;
	}
	//std::cout<< " Possible external species in curly bracket and internal species in () is shown below \n";
	if (!possible_species_1.empty()) {
		for (int i = 0; i<possible_species_1.size();i++){

			mband::assign_label(graph,int_vector[0],possible_species_1[i]);

			std::vector<int> initial_species_2 = mband::generate_edge_species(graph, int_vector[1]);

			std::vector<std::vector<int>> possible_species_2 = mband::findmatch(initial_species_2,int_vector[1]);
			//end_code
			if (!possible_species_2.empty()){
				for (int j = 0; j < possible_species_2.size();j++){
				mband::assign_label(graph,int_vector[1],initial_species_2);				
				mband::assign_label(graph,int_vector[1],possible_species_2[j]);
				interaction_species.push_back({possible_species_1[i],possible_species_2[j]});
				ext_legs.push_back(mband::external_species(graph));
				
				std::cout<< fermionic_edge.size();
				std::vector<int> v;
				//std::cout<<"(";
				for (int x =0; x < fermionic_edge.size(); x++){
					v.push_back(graph[fermionic_edge[x]].g_struct_.species_);
					//std::cout<<graph[fermionic_edge[x]].g_struct_.species_;				
				}
				fermionic_species.push_back(v);
				v.clear();
				//std::cout<<")"<<std::endl;			    
				}
				mband::assign_label(graph,int_vector[1],initial_species_2);
				
			}
			else {mband::assign_label(graph,int_vector[1],initial_species_2);}
			
		mband::reset_species(graph,int_vector);	
		}
	}
	else {std::cout <<"No match found" <<std::endl;}

mband::print_assigned_species(interaction_species);

}


void mband::solve_multiband_1(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs){
    AmiGraph::edge_vector_t bvector;
	std::vector<AmiGraph::edge_vector_t> int_vector;
	g.find_bosonic_edges(graph,bvector);
	g.find_internal_fermionic_edges(graph,fermionic_edge);
	for (int i=0;i <bvector.size();i++){
		//g.print_edge_info(bvector[i],gself2);
		std::cout << "Edge = ("  << graph[source(bvector[i],graph)].index_<<"," << graph[target(bvector[i],graph)].index_ <<")" 
		<<std::endl;
	}
	for (auto b: bvector){
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
	}	
	mband::find_interaction(graph,bvector,int_vector);
	mband::print_interactions(graph,bvector, int_vector);
	std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);
	std::vector<std::vector<int>> possible_species_1;
	if (mband::Hartee_fock){
		possible_species_1 = findmatch(initial_species_1,int_vector[0]);
	}
	else{
		possible_species_1 = mband::interaction_legs;
	}
	std::cout<< " Possible external species in curly bracket and internal species in () is shown below \n";
	if (!possible_species_1.empty()) {
		for (int i = 0; i<possible_species_1.size();i++){

		mband::assign_label(graph,int_vector[0],possible_species_1[i]);
        interaction_species.push_back({possible_species_1[i]});
		ext_legs.push_back(mband::external_species(graph));
		std::vector<int> v;
		std::cout<<"(";
		for (int x =0; x < fermionic_edge.size(); x++){
			v.push_back(graph[fermionic_edge[x]].g_struct_.species_);
			std::cout<<graph[fermionic_edge[x]].g_struct_.species_;				
		}
		fermionic_species.push_back(v);
		v.clear();
		std::cout<<")"<<std::endl;			    
		}
		
		mband::reset_species(graph,int_vector);	
		
	}
	else {std::cout <<"No match found" <<std::endl;}

mband::print_assigned_species(interaction_species);

  
}





void  mband::generate_eps_alpha(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge, std::vector<AmiBase::epsilon_t> &ept, 
std::vector<AmiBase::alpha_t>  &alpha) {
	
	for (int i =0; i < fermionic_edge.size();i++){
		ept.push_back(graph[fermionic_edge[i]].g_struct_.eps_);
		alpha.push_back(graph[fermionic_edge[i]].g_struct_.alpha_);
		
		
	}

}




void mband::solve_multiband_5(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs){
    AmiGraph::edge_vector_t bvector;
	std::vector<AmiGraph::edge_vector_t> int_vector;
	g.find_bosonic_edges(graph,bvector);
	g.find_internal_fermionic_edges(graph,fermionic_edge);
	for (int i=0;i <bvector.size();i++){
		//g.print_edge_info(bvector[i],gself2);
		std::cout << "Edge = ("  << graph[source(bvector[i],graph)].index_<<"," << graph[target(bvector[i],graph)].index_ <<")" 
		<<std::endl;
	}
	for (auto b: bvector){
		bosonic_Alpha.push_back(graph[b].g_struct_.alpha_);
	}
	mband::find_interaction(graph,bvector,int_vector);
	mband::print_interactions(graph,bvector, int_vector);
	std::vector<std::vector<int>> possible_species_1;
	std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);
	if (mband::Hartee_fock){
		possible_species_1 = findmatch(initial_species_1,int_vector[0]);
	}
	else{
		possible_species_1 = mband::interaction_legs;
	}
	
	std::cout<< " Possible external species in curly bracket and internal species in () is shown below \n";
	if (!possible_species_1.empty()){
		for (int i = 0; i<possible_species_1.size();i++){
			
			std::vector<int> initial_species_1 = mband::generate_edge_species(graph, int_vector[0]);	
			mband::assign_label(graph,int_vector[0],possible_species_1[i]);
			std::vector<int> initial_species_2 = mband::generate_edge_species(graph, int_vector[1]);
			std::vector<std::vector<int>> possible_species_2 = mband::findmatch(initial_species_2,int_vector[1]);

			if (!possible_species_2.empty()){
				for (int j = 0; j < possible_species_2.size();j++){
				//mband::print_match(possible_species_2[j],2);
				mband::assign_label(graph,int_vector[1],initial_species_2);				
				mband::assign_label(graph,int_vector[1],possible_species_2[j]);
				std::vector<int> initial_species_3 = mband::generate_edge_species(graph, int_vector[2]);
				std::vector<std::vector<int>> possible_species_3 = mband::findmatch(initial_species_3,int_vector[2]);
				
					if (!possible_species_3.empty()){
						for (int k = 0; k < possible_species_3.size();k++){
							//mband::print_match(possible_species_3[k],3);
							mband::assign_label(graph,int_vector[2],initial_species_3);
							mband::assign_label(graph,int_vector[2],possible_species_3[k]);											
							std::vector<int> initial_species_4 = mband::generate_edge_species(graph, int_vector[3]);
							std::vector<std::vector<int>> possible_species_4 = mband::findmatch(initial_species_4,int_vector[3]);					
							if (!possible_species_4.empty()){
								for (int l = 0; l < possible_species_4.size();l++){
									//num++
								//mband::print_match(possible_species_4[l],4);
								mband::assign_label(graph,int_vector[3],initial_species_4);							
								mband::assign_label(graph,int_vector[3],possible_species_4[l]);
								std::vector<int> initial_species_5 = mband::generate_edge_species(graph, int_vector[4]);
								std::vector<std::vector<int>> possible_species_5 = mband::findmatch(initial_species_5,int_vector[4]);
										if (!possible_species_5.empty()){
											for (int m = 0; m < possible_species_5.size();m++){
												mband::assign_label(graph,int_vector[4],initial_species_5);							
												mband::assign_label(graph,int_vector[4],possible_species_5[m]);
											    interaction_species.push_back({possible_species_1[i],possible_species_2[j],possible_species_3[k],possible_species_4[l],possible_species_5[m]});
												ext_legs.push_back(mband::external_species(graph));
												std::cout<< fermionic_edge.size();
												
													std::cout<<"(";		
												std::vector<int> v;
												for (int x =0; x < fermionic_edge.size(); x++){
													v.push_back(graph[fermionic_edge[x]].g_struct_.species_);
													std::cout<<graph[fermionic_edge[x]].g_struct_.species_;
												}
												fermionic_species.push_back(v);
												v.clear();
												std::cout<<")"<<std::endl;	
											}
										}
								else { mband::assign_label(graph,int_vector[4],initial_species_5);}				
							}mband::assign_label(graph,int_vector[3],initial_species_4);
								
							}
							else { mband::assign_label(graph,int_vector[3],initial_species_4);}				
						}mband::assign_label(graph,int_vector[2],initial_species_3);
										
					}
					else{mband::assign_label(graph,int_vector[2],initial_species_3);}
				}mband::assign_label(graph,int_vector[1],initial_species_2);			
			}
			else {mband::assign_label(graph,int_vector[1],initial_species_2);}		
		mband::reset_species(graph,int_vector);	}
	}
	else{std::cout<< "NO match found \n";}
mband::print_assigned_species(interaction_species);
}




	
void mband::solve_multiband(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t &fermionic_edge,std::vector<std::vector<int>> &fermionic_species,std::vector<std::vector<std::vector<int>>> &interaction_species,std::vector<std::vector<int>> &bosonic_Alpha,std::vector<std::vector<int>> &ext_legs){
int ord = g.graph_order(graph);
if (ord == (int) 1){
 mband::solve_multiband_1(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,ext_legs);}
else if (ord == (int) 2){
 mband::solve_multiband_2(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,ext_legs);}
else if (ord == (int) 3){
 mband::solve_multiband_3(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,ext_legs);}
else if (ord ==(int) 4){
 mband::solve_multiband_4(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,ext_legs);}
 else if (ord ==(int) 5){
 mband::solve_multiband_5(graph,fermionic_edge,fermionic_species,interaction_species,bosonic_Alpha,ext_legs);}
else{
std::cout << "this order not possible";}
}





		


		                 	






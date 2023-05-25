#include "mini_ami.hpp"
#include "algorithm"


void mband::print_assigned_species(std::vector<std::vector<std::vector<int>>> interaction_species){
	     for (int i = 0; i < interaction_species.size(); i++) {
        std::cout << "Interactions are  " << i << ":" << std::endl;
        // Iterate through each row
        for (int j = 0; j < interaction_species[i].size(); j++) {
            // Iterate through each column
            for (int k = 0; k < interaction_species[i][j].size(); k++) {
                std::cout << interaction_species[i][j][k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
	
	}
void mband::print_interactions(AmiGraph::graph_t &graph,AmiGraph::edge_vector_t b_vector, std::vector<AmiGraph::edge_vector_t> f_vector){
	std::cout  << " Printing interaction edges of each bosonic line " << std::endl;
	std::cout << " (1) incoming source, (2) outgoing source , (3)  incoming target , (4) outgoing target" << "\n ";
    //std::assert(b_vector.size()==f_vector.size());
	 for (int i=0;i <b_vector.size();i++){
		 std::cout << "Edge = ("  << graph[source(b_vector[i],graph)].index_<<"," << graph[target(b_vector[i],graph)].index_ <<") \n "  ;
		 
		 for(int j =0; j<f_vector[i].size();j++)
		 {
		 std::cout << "("<<j+1 << ")" << "edge is " <<  "("<<graph[source(f_vector[i][j],graph)].index_ << ","<<graph[target(f_vector[i][j],graph)].index_  
		 << ")" <<std::endl;
			 
		
	 
	
	     }
	 }
}
void mband::print_match(std::vector<int> vec,int ord) {
	if (vec.empty()){
		std::cout<<ord<<"empty" << std::endl;
	}
	else{
	std::cout<<"for order " << ord  << std::endl;
    for (int elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
	}
}	
	

double mband::perturb(double number) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.01, 0.9);
    std::uniform_int_distribution<int> sign_dist(0, 1);
    int sign = (sign_dist(gen) == 0) ? -1 : 1;
    return number + (dist(gen) * 1e-5 * sign);
}


std::vector<double> mband::sumVectors(std::vector<std::vector<double>> vectors)
{
    std::vector<double> sum(vectors[0].size(), 0);
    for (auto vec : vectors)
    {
        for (int i = 0; i < vec.size(); i++)
        {   
			if (sum[i] == 0){
            sum[i] += vec[i];
			}
        }
    }
    return sum;
}
std::vector<std::complex<double>> mband::convertToComplex(const std::vector<double> vec) {
    std::vector<std::complex<double>> cplx_vec;
    cplx_vec.reserve(vec.size()); 

    for (auto elem : vec) {
        cplx_vec.emplace_back(elem, 0.0);
    }

    return cplx_vec;
}

std::vector<std::vector<double>>  mband::band_to_hab(std::vector<std::vector<int>> band) {
	std::vector<double> hab = mband::energy;
    std::vector<std::vector<double>> band_hab(band.size(), std::vector<double>(band[0].size()));
    for (int i = 0; i < band.size(); i++) {
        for (int j = 0; j < band[i].size(); j++) {
            band_hab[i][j] =  -1*mband::perturb(hab[band[i][j] - 1]);
			//band_hab[i][j] =  -1*hab[band[i][j] - 1];
        }
    }
    std::cout << "Band index replaced by energy:\n";
    for (int i = 0; i < band_hab.size(); i++)
    {
        for (int j = 0; j < band_hab[i].size(); j++)
        {
            std::cout << band_hab[i][j] << " ";
        }
        std::cout << "\n\n";
    }
    return band_hab;
}

std::vector<std::complex<double>> mband::generate_ept(std::vector<std::vector<int>> epsilon, std::vector<double> band_value) {
	//std::sort(epsilon.begin(), epsilon.end());
    //epsilon.erase(std::unique(epsilon.begin(), epsilon.end()), epsilon.end());
    std::vector<std::vector<double>> results;
    std::vector<double> collect;
    for (int  i = 0; i < epsilon.size(); i++) {
        
        for (int j = 0; j < epsilon[i].size(); j++) {
            collect.push_back(epsilon[i][j] * band_value[i]);
        }
        results.push_back(collect);
        collect.clear();;

    }
    std::vector<double> sum = mband::sumVectors(results);



    return mband::convertToComplex(sum);
    
}


double mband::Umatch(const std::vector<std::vector<int>>& int_matrix, const std::vector<double>& int_value, const std::vector<std::vector<int>>& int_species) {
    double U = 1.00;
    for (int i = 0; i < int_species.size(); i++) {
        for (int j = 0; j < int_matrix.size(); j++) {
            if ((int_species[i][1] == int_matrix[j][1] && int_species[i][2] == int_matrix[j][2] && int_species[i][3]
                == int_matrix[j][3] && int_species[i][0] == int_matrix[j][0])) {
                U = U * int_value[j];
            }
        }
       }
    return U;
}
/*
double mband::Umatch(const std::vector<std::vector<int>>& int_matrix, const std::vector<double>& int_value, 
                     const std::vector<std::vector<int>>& int_species) {
    double U = 1.00;
    std::unordered_map<int, double> int_map;
    
    // Create a map to store the interaction matrix values
    for (int i = 0; i < int_matrix.size(); i++) {
        int key = int_matrix[i][0] * 1000 + int_matrix[i][1] * 100 + int_matrix[i][2] * 10 + int_matrix[i][3];
        int_map[key] = int_value[i];
    }
    
    // Iterate over the input vector and lookup the interaction values in the map
    for (const auto& species : int_species) {
        int key = species[0] * 1000 + species[1] * 100 + species[2] * 10 + species[3];
        auto it = int_map.find(key);
        if (it != int_map.end()) {
            U *= it->second;
            continue; // Terminate the inner loop early
        }
    }
    
    return U;
}
*/

void mband::filter(std::vector<std::vector<int>>& possible_species, const std::vector<int>& list) {
    if (list.empty()) {

    }
    else {
        size_t x = static_cast<size_t>(list[0]);
        size_t y = static_cast<size_t>(list[1]);

        std::vector<std::vector<int>> filter;
        filter.reserve(possible_species.size());

        for (const auto& vec : possible_species) {
            if (vec[x] == vec[y]) {
                filter.emplace_back(vec);
            }
        }

        possible_species = std::move(filter);
    }

}

std::vector<int>  mband::interaction_index(const  std::vector<std::vector<int>>& int_species) {
   std::vector<int> vec;
   const std::vector<std::vector<int>> int_matrix = mband::interaction_legs;
    for (int i = 0; i < int_species.size(); i++) {
        for (int j = 0; j < int_matrix.size(); j++) {
            if ((int_species[i][1] == int_matrix[j][1] && int_species[i][2] == int_matrix[j][2] && int_species[i][3]
                == int_matrix[j][3] && int_species[i][0] == int_matrix[j][0])) {
                vec.push_back(j);
            }
        }
        
    }
   
    return vec;
   
}

/*
std::vector<int> mband::interaction_index(const std::vector<std::vector<int>>& int_species) {
    std::vector<int> vec;
    const std::vector<std::vector<int>>& int_matrix = mband::interaction_legs;
    std::unordered_map<int, std::vector<int>> int_map;
    
    // Create a map to store the interaction matrix
    for (int i = 0; i < int_matrix.size(); i++) {
        int key = int_matrix[i][0] * 1000 + int_matrix[i][1] * 100 + int_matrix[i][2] * 10 + int_matrix[i][3];
        int_map[key].push_back(i);
    }
    
    // Iterate over the input vector and lookup the interaction indices in the map
    for (const auto& species : int_species) {
        int key = species[0] * 1000 + species[1] * 100 + species[2] * 10 + species[3];
        auto it = int_map.find(key);
        if (it != int_map.end()) {
            vec.insert(vec.end(), it->second.begin(), it->second.end());
        }
    }
    
    return vec;
}

*/





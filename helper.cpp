#include "mini_ami.hpp"

double mband::perturb(double number) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.01, 0.9);
    std::uniform_int_distribution<int> sign_dist(0, 1);
    int sign = (sign_dist(gen) == 0) ? -1 : 1;
    return number + (dist(gen) * 1e-2 * sign/2);
}


std::vector<double> mband::sumVectors(std::vector<std::vector<double>> vectors)
{
    std::vector<double> sum(vectors[0].size(), 0);
    for (auto vec : vectors)
    {
        for (int i = 0; i < vec.size(); i++)
        {
            sum[i] += vec[i];
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
            //band_hab[i][j] =  mband::perturb(hab[band[i][j] - 1]);
			band_hab[i][j] =  hab[band[i][j] - 1];
        }
    }
    std::cout << "Band index replaced by energy:\n";
    for (int i = 0; i < band_hab.size(); i++)
    {
        for (int j = 0; j < band_hab[i].size(); j++)
        {
            std::cout << band_hab[i][j] << " ";
        }
        std::cout << "\n";
    }
    return band_hab;
}

std::vector<std::complex<double>> mband::generate_ept(std::vector<std::vector<int>> epsilon, std::vector<double> band_value) {
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


double mband::Umatch(std::vector<std::vector<int>> int_matrix, std::vector<double> int_value, std::vector<std::vector<int>> int_species) {
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

void filter(std::vector<std::vector<int>>& possible_species, const std::vector<double>& list) {
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

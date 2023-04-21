#include "mini_ami.hpp"

void mband::changeEqualNumbers(std::vector<double>& vec) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1, 1);
    // Iterate over the vector
    for (int i = 0; i < vec.size(); i++) {
        // Check if the current element is equal to the next element
        if (i + 1 < vec.size() && vec[i] == vec[i + 1]) {
            // Find the end of the group of equal numbers
            int j = i + 1;
            while (j < vec.size() && vec[j] == vec[i]) {
                j++;
            }
            j--;

            // Reduce the values of the group by 0.1
            for (int k = i; k <= j; k++) {
                vec[k] = vec[k]+ dist(gen)*1e-7;
            }

            // Move to the next element after the group
            i = j;
        }
    }
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
	changeEqualNumbers(sum);


    return mband::convertToComplex(sum);
    
}
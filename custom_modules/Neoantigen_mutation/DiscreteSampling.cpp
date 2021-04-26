#include <iostream>
#include <random>
#include <vector>
#include <array>

int main()
{
    const int sampleSize = 3;   // Size of the sample
    std::vector<double> weights = {47,38,15}; // 10 possible outcome with different weights

    std::random_device rd;
    std::mt19937 generator(rd());

    /// WITH REPLACEMENT

    std::discrete_distribution<int> distribution(weights.begin(), weights.end());

    std::array<int, sampleSize> p ={};
    for(int i=0; i< 200; ++i){
        int number = distribution(generator);
        ++p[number];
    }

    std::cout << "Discrete_distribution with replacement:" << std::endl;
    for (int i=0; i < 3; ++i)
    std::cout << i << ": " << std::string(p[i],'*') << std::endl;

    for (int i=0; i < 3; ++i)
    std::cout << i << ": " << p[i] << std::endl;

    /// WITHOUT REPLACEMENT
    //
    // p = {};
    // for(int i=0; i<sampleSize; ++i){
    //     std::discrete_distribution<int> distribution(weights.begin(), weights.end());
    //     int number = distribution(generator);
    //     weights[number] = 0; // the weight associate to the sampled value is set to 0
    //     ++p[number];
    // }
    //
    // std::cout << "Discrete_distribution without replacement:" << std::endl;
    // for (int i=0; i<10; ++i)
    // std::cout << i << ": " << std::string(p[i],'*') << std::endl;


    return 0;
}

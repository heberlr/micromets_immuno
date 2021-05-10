#include <iostream>
#include <random>
#include <vector>
#include <array>

class Cell{
  public:
    int ID;
    std::vector<int> neoantigen_signature;
};

void print(std::vector<int> const input)
{
  for (int i = 0; i < input.size(); i++)
    std::cout << input[i] << ' ';
}

int main()
{
    const int sampleSize = 3;   // Size of the sample
    std::vector<double> weights = {47,38,15}; // 10 possible outcome with different weights

    std::vector<int> neoant_sig(10,0.0);
    std::vector<Cell*> all_cells;
    Cell pCell;
    pCell.ID = 0;
    pCell.neoantigen_signature = neoant_sig;
    all_cells.push_back(&pCell);
    std::cout << "Size: " << all_cells.size() << " ID: " << pCell.ID << " Neoatigen sig: ";
    print(pCell.neoantigen_signature);
    std::cout << std::endl;

    // std::random_device rd;
    // std::mt19937 generator(rd());
    //
    // /// WITH REPLACEMENT
    //
    // std::discrete_distribution<int> distribution(weights.begin(), weights.end());
    //
    // std::array<int, sampleSize> p ={};
    // for(int i=0; i< 200; ++i){
    //     int number = distribution(generator);
    //     ++p[number];
    // }
    //
    // std::cout << "Discrete_distribution with replacement:" << std::endl;
    // for (int i=0; i < 3; ++i)
    // std::cout << i << ": " << std::string(p[i],'*') << std::endl;
    //
    // for (int i=0; i < 3; ++i)
    // std::cout << i << ": " << p[i] << std::endl;

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

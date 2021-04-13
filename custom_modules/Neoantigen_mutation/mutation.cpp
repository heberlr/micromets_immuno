#include <iostream>
#include <vector>
#include <random>
#include <fstream>

std::random_device rd;
std::mt19937 gen(rd());

double UniformRandom(){
	return std::generate_canonical<double, 10>(gen);
}

void mutation (std::vector<bool> &Sequence, double prob_mutation ){
  for (int i=0; i < Sequence.size(); i++)
    if (UniformRandom() < prob_mutation) Sequence[i] = !Sequence[i];
}

void printVector (const std::vector<bool> &Sequence){
  std::cout <<"[ ";
  for (int i=0; i < Sequence.size(); i++)
    std::cout << Sequence[i] << " ";
  std::cout <<"]"<< std::endl;
}

int hammingDist(const std::vector<bool> &Seq1, const std::vector<bool> &Seq2){
  if (Seq1.size() != Seq1.size()){ std::cout << "Error: size of Seq1 uncompatible to Seq2" << std::endl; exit(0);}
  int count = 0;
  for (int i=0; i < Seq1.size(); i++){
    if (Seq1[i] != Seq2[i]) count++;
  }
  return count;
}



int main( int argc, char* argv[] )
{
  if( argc != 3 ) std::cout << "Error: give 2 arguments: SizeSeq and file" << std::endl;
  const double Size = atof(argv[1]);
  std::ofstream File;
  File.open(argv[2]);
  std::vector<int> Generation = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  std::vector<double> MutationRate = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  std::vector<bool> BooleanVectorRef(Size,false);
  std::vector<bool> BooleanVector = BooleanVectorRef;

  for (int i=0; i < MutationRate.size(); i++){
    for (int j=0; j < Generation.size(); j++){
      mutation(BooleanVector,MutationRate[i]);
      File << Size << " " << MutationRate[i] << " " << Generation[j] << " " << hammingDist(BooleanVectorRef,BooleanVector) << std::endl;
    }
    BooleanVector = BooleanVectorRef;
  }
  File.close();

  return 0;
}

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

#include "./submodel_data_structures.h"

#ifndef __external_immune__
#define __external_immune__

extern Submodel_Information external_immune_info;

void external_immune_model_setup( void );

void external_immune_model( double dt );

class antigen{
	public:
		int ID;
		std::vector<double> signature;
		double time;
		int probability_weight;
};

class AntigenLibrary{
	public:
			std::vector<antigen> Collection;
			std::vector<antigen> TempCollection;

			//Methods
			bool check_library(const std::vector<double> &);
			void add_antigen(const std::vector<double> &, const double);
			void update_collection(const double, const double);
			std::vector<double> sample_antigen();
			void print();
};

#endif

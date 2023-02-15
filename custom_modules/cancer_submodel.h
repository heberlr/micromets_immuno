#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

#include "./submodel_data_structures.h"

#include "./immune_submodels.h"

#ifndef __cancer_submodel__
#define __cancer_submodel__

extern Submodel_Information cancer_submodel_info;

void cancer_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt );
void cancer_phenotype( Cell* pCell, Phenotype& phenotype, double dt );
void cancer_mechanics( Cell* pCell, Phenotype& phenotype, double dt );
void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt );
void cancer_submodel_setup( void );

#endif

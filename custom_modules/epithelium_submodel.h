#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

#include "./submodel_data_structures.h"

#ifndef __epithelium_submodel__
#define __epithelium_submodel__

extern Submodel_Information epithelium_submodel_info;

void epithelium_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt );
void epithelium_phenotype( Cell* pCell, Phenotype& phenotype, double dt );
void epithelium_mechanics( Cell* pCell, Phenotype& phenotype, double dt );
void epithelium_submodel_setup( void );

#endif

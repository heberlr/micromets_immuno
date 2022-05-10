#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

#include "./submodel_data_structures.h"

#include "./external_immune.h"

#ifndef __immune_submodels__
#define __immune_submodels__

/* functions for checking nearby cells */
Cell* check_for_live_neighbor_for_interaction( Cell* pAttacker , double dt );
Cell* check_for_dead_neighbor_for_interaction( Cell* pAttacker , double dt );

/* functions for cell-cell adhesion */
bool attempt_immune_cell_attachment( Cell* pAttacker, Cell* pTarget , double dt );
Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt );

/* functions for vascularization */
std::vector<double> choose_vascularized_position( void );
void choose_initialized_voxels( void );

/* functions for immune cells */
void CD8_Tcell_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2 , double dt );
void CD8_Tcell_phenotype( Cell* pCell, Phenotype& phenotype, double dt );
void CD8_Tcell_mechanics( Cell* pCell, Phenotype& phenotype, double dt );
void macrophage_phenotype( Cell* pCell, Phenotype& phenotype, double dt );
void macrophage_mechanics( Cell* pCell, Phenotype& phenotype, double dt );
void DC_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2 , double dt );
void DC_phenotype( Cell* pCell, Phenotype& phenotype, double dt );
void DC_mechanics( Cell* pCell, Phenotype& phenotype, double dt );
void CD4_Tcell_mechanics( Cell* pCell, Phenotype& phenotype, double dt );
void CD4_Tcell_phenotype( Cell* pCell, Phenotype& phenotype, double dt );

/* chemotaxis of immune cells */
void immune_cell_motility_direction( Cell* pCell, Phenotype& phenotype , double dt );

/* check for cells out of bound */
bool check_for_out_of_bounds( Cell* pC , double tolerance );

/* create cells at the beginning */
void initial_immune_cell_placement( void );
void create_infiltrating_immune_cell_initial( Cell_Definition* pCD );

// immune cell recruitment
void immune_cell_recruitment( double dt );
void create_infiltrating_immune_cell( Cell_Definition* pCD );
void create_infiltrating_immune_cell( std::string cell_name );
void create_infiltrating_Tcell( void );
void create_infiltrating_macrophage( void );
void create_infiltrating_DCcell( void );
void create_infiltrating_CD4Tcell( void );

// register models
void immune_submodels_setup( void );

int Hamming_Distance(const std::vector<double> &Seq1, const std::vector<double> &Seq2);

#endif

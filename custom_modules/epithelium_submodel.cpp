#include "./epithelium_submodel.h"

using namespace PhysiCell;

std::string epithelium_submodel_version = "0.4.0";

Submodel_Information epithelium_submodel_info;

void epithelium_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt )
{
	// elastic adhesions
	standard_elastic_contact_function( pC1,p1, pC2, p2, dt );

	return;
}

void epithelium_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	static int proinflammatory_cytokine_index = microenvironment.find_density_index("pro-inflammatory cytokine");
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" );

	phenotype.motility.is_motile = false;

	// if I am dead, remove all adhesions
	if( phenotype.death.dead == true )
	{
		// detach all attached cells
		// remove_all_adhesions( pCell );

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
	}

	int cycle_G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int cycle_S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
	// turn off proliferation
	pCell->phenotype.cycle.data.transition_rate(cycle_G0G1_index,cycle_S_index) = 0.0;

	// Normal cells neighbor to melanoma cells dead by inflammation (kill cells with high pro-inf concentration)
	// if( (pCell->nearest_density_vector())[ proinflammatory_cytokine_index ] > 0.9 ) // Check neighbor if have high density of melanoma cells.
	// {
	// 		pCell->start_death( apoptosis_index );
	// }

	// if I am dead, don't bother executing this function again
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
	}

	return;
}

void epithelium_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");

	//melanoma cells are movable
	// if ( pCell->custom_data["melanoma_cell_chemokine_secretion_activated"] > 0.1)
	// 	pCell->is_movable = true;
	// else
	// 	pCell->is_movable = false;

	// if I'm dead, don't bother
	if( phenotype.death.dead == true )
	{
		// the cell death functions don't automatically turn off custom functions,
		// since those are part of mechanics.
		// remove_all_adhesions( pCell );

		// Let's just fully disable now.
		pCell->functions.custom_cell_rule = NULL;
		pCell->functions.contact_function = NULL;

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
		return;
	}

	// this is now part of contact_function
	/*
	// if I'm adhered to something ...
	if( pCell->state.neighbors.size() > 0 )
	{
	// add the elastic forces
	extra_elastic_attachment_mechanics( pCell, phenotype, dt );
}
*/
return;
}

void epithelium_submodel_setup( void )
{
	Cell_Definition* pCD;

	// set up any submodels you need

	// receptor trafficking
	//receptor_dynamics_model_setup(); // done
	// pathogen replication
	//internal_pathogen_model_setup();
	// single-cell response
	// internal_pathogen_response_model_setup();

	// set up epithelial cells
	// set version info
	epithelium_submodel_info.name = "epithelium model";
	epithelium_submodel_info.version = epithelium_submodel_version;
	// set functions
	epithelium_submodel_info.main_function = NULL;
	epithelium_submodel_info.phenotype_function = epithelium_phenotype;
	epithelium_submodel_info.mechanics_function = epithelium_mechanics;

	// what microenvironment variables do you expect?
	epithelium_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" );
	epithelium_submodel_info.microenvironment_variables.push_back( "chemokine" );

	// what custom data do I need?
	//epithelium_submodel_info.cell_variables.push_back( "something" );
	// register the submodel
	epithelium_submodel_info.register_model();
	// set functions for the corresponding cell definition
	pCD = find_cell_definition( "lung cell" );
	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
	pCD->functions.contact_function = epithelium_contact_function;

	return;
}

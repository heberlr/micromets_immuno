#include "./melanoma_submodel.h"

using namespace PhysiCell;

std::string melanoma_submodel_version = "0.1.0";

Submodel_Information melanoma_submodel_info;

// Cells that will receive a nudge to come back to the domain
extern std::vector<Cell*> cells_to_move_from_edge;

void melanoma_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt )
{
	// elastic adhesions
	standard_elastic_contact_function( pC1,p1, pC2, p2, dt );

	return;
}

double strain_based_proliferation( Cell* pCell )
{
	static double  max_pressure = 10.0; // maximum tolerated pressure of melanoma cell (proliferation) [10.0 microns]
	if( pCell->state.simple_pressure < max_pressure )
	{
		return pow( (max_pressure - pCell->state.simple_pressure)/max_pressure, 1.0 );
	}
	return 0.0;
}

void melanoma_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" );

	phenotype.motility.is_motile = false;

	// T-cell based death
	TCell_induced_apoptosis(pCell, phenotype, dt );

	// if I am dead, remove all adhesions
	if( phenotype.death.dead == true )
	{
		// detach all attached cells
		// remove_all_adhesions( pCell );
		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
	}

	// Mechanical contribution to proliferation
	double mechanics_factor = strain_based_proliferation( pCell );

	// proliferation rate based on mechanical aspect
	int cycle_G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int cycle_S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
	pCell->phenotype.cycle.data.transition_rate(cycle_G0G1_index,cycle_S_index) = mechanics_factor*pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(cycle_G0G1_index,cycle_S_index);

	// if I am dead, don't bother executing this function again
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
	}

	return;
}

void melanoma_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");

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

	// bounds check
	if( check_for_out_of_bounds( pCell , 10.0 ) )
	{
		#pragma omp critical
		{ cells_to_move_from_edge.push_back( pCell ); }
	}

	return;
}

void melanoma_submodel_setup( void )
{
	Cell_Definition* pCD;

	// set up any submodels you need

	// set up epithelial cells
	// set version info
	melanoma_submodel_info.name = "melanoma model";
	melanoma_submodel_info.version = melanoma_submodel_version;
	// set functions
	melanoma_submodel_info.main_function = NULL;
	melanoma_submodel_info.phenotype_function = melanoma_phenotype;
	melanoma_submodel_info.mechanics_function = melanoma_mechanics;

	// what microenvironment variables do you expect?
	melanoma_submodel_info.microenvironment_variables.push_back( "debris" );

	// what custom data do I need?
	//melanoma_submodel_info.cell_variables.push_back( "something" );
	// register the submodel
	melanoma_submodel_info.register_model();
	// set functions for the corresponding cell definition
	pCD = find_cell_definition( "melanoma cell" );
	pCD->functions.update_phenotype = melanoma_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = melanoma_submodel_info.mechanics_function;
	pCD->functions.contact_function = melanoma_contact_function;

	return;
}

void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" );
	static int debris_index = microenvironment.find_density_index( "debris" );

	if( pCell->custom_data["TCell_contact_time"] > pCell->custom_data["TCell_contact_death_threshold"] )
	{
		#pragma omp critical
		{
			std::cout << "\t\t\t\t" << pCell << " (of type " << pCell->type_name <<  ") died from T cell contact" << std::endl;
		}

		// induce death
		pCell->start_death( apoptosis_index );
		pCell->phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];

		pCell->functions.update_phenotype = NULL;
	}

	return;
}

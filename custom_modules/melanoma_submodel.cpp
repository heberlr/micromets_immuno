#include "./melanoma_submodel.h"

using namespace PhysiCell;

std::string melanoma_submodel_version = "0.1.0";

Submodel_Information melanoma_submodel_info;

void melanoma_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2, double dt )
{
	// elastic adhesions
	standard_elastic_contact_function( pC1,p1, pC2, p2, dt );

	return;
}

void melanoma_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	static int proinflammatory_cytokine_index = microenvironment.find_density_index("pro-inflammatory cytokine");
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

	int cycle_G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int cycle_S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
	//proliferation, mechanics, and chemokine secretion activated

	// Mechanical contribution to proliferation
	double mechanics_factor = pow(500.0 - pCell->state.simple_pressure,0.5)/pow(500.0,0.5);
	if (mechanics_factor < 0.0) mechanics_factor = 0.0;
	// proliferation rate based on mechanical aspect
	pCell->phenotype.cycle.data.transition_rate(cycle_G0G1_index,cycle_S_index) = parameters.doubles("prolif_rate_CancerCell")*mechanics_factor;
	pCell->custom_data["melanoma_cell_chemokine_secretion_activated"] = 1.0;

	static int chemokine_index = microenvironment.find_density_index( "chemokine" );
	if( pCell->custom_data["melanoma_cell_chemokine_secretion_activated"] > 0.1 && phenotype.death.dead == false )
	{
		double rate = 1.0; //AV; // P;
		rate /= pCell->custom_data["max_apoptosis_half_max"];
		if( rate > 1.0 )
		{ rate = 1.0; }
		rate *= pCell->custom_data[ "melanoma_cell_chemokine_secretion_rate" ];

		phenotype.secretion.secretion_rates[chemokine_index] = rate;
		phenotype.secretion.saturation_densities[chemokine_index] = 1.0;
	}

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

void melanoma_submodel_setup( void )
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
	melanoma_submodel_info.name = "melanoma model";
	melanoma_submodel_info.version = melanoma_submodel_version;
	// set functions
	melanoma_submodel_info.main_function = NULL;
	melanoma_submodel_info.phenotype_function = melanoma_phenotype;
	melanoma_submodel_info.mechanics_function = melanoma_mechanics;

	// what microenvironment variables do you expect?
	melanoma_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" );
	melanoma_submodel_info.microenvironment_variables.push_back( "chemokine" );
	melanoma_submodel_info.microenvironment_variables.push_back( "anti-inflammatory cytokine" );

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
	static int proinflammatory_cytokine_index = microenvironment.find_density_index("pro-inflammatory cytokine");
	static int antiinflammatory_cytokine_index = microenvironment.find_density_index("anti-inflammatory cytokine");

	if( pCell->custom_data["TCell_contact_time"] > pCell->custom_data["TCell_contact_death_threshold"] )
	{
		// make sure to get rid of all adhesions!
		// detach all attached cells
		// remove_all_adhesions( pCell );

		#pragma omp critical
		{
			std::cout << "\t\t\t\t" << pCell << " (of type " << pCell->type_name <<  ") died from T cell contact" << std::endl;
		}

		// induce death
		pCell->start_death( apoptosis_index );

		pCell->phenotype.secretion.secretion_rates[proinflammatory_cytokine_index] = 0;
		pCell->phenotype.secretion.secretion_rates[antiinflammatory_cytokine_index] = pCell->custom_data["antiinflammatory_cytokine_secretion_rate"];
		pCell->phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];

		pCell->functions.update_phenotype = NULL;
	}

	return;
}

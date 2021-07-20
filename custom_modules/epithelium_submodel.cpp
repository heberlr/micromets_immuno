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

bool strain_based_apoptosis( Cell* pCell )
{
	static int strain_index = pCell->custom_data.find_variable_index( "mechanical_strain" );
	static double  max_strain = 0.75; //<max_mechanical_strain description="maximum tolerated deformation of lung cell (death)" type="double" units="micron">0.75</max_mechanical_strain>
	if( pCell->custom_data[strain_index] <= max_strain ) return false;

	std::vector<Cell*> neighbors = pCell->cells_in_my_container();//find cells in a neighbourhood of melanoma cells
	static int melanoma_cell_type = get_cell_definition( "melanoma cell" ).type;
	int n = 0;
	Cell* pTestCell;
	while( n < neighbors.size() )
	{
		pTestCell = neighbors[n];
		if ( pTestCell != pCell && pTestCell->phenotype.death.dead == false && pTestCell->type == melanoma_cell_type){
			double cell_cell_distance = sqrt((pTestCell->position[0]-pCell->position[0])*(pTestCell->position[0]-pCell->position[0])+(pTestCell->position[1]-pCell->position[1])*(pTestCell->position[1]-pCell->position[1]));
			double radius_DC = pCell->phenotype.geometry.radius; // (Adrianne) radius of DC)
			double radius_test_cell = pTestCell->phenotype.geometry.radius;
			if( cell_cell_distance <= parameters.doubles("epsilon_distance")*(radius_DC+radius_test_cell) ){
				return true;
			}
		}
		n++;
	}
	return false;
}

void epithelium_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" );

	phenotype.motility.is_motile = false;

	// update mechanical strain
	static int strain_index = pCell->custom_data.find_variable_index( "mechanical_strain" );
	static int ECM_attachment_point_index = pCell->custom_data.find_vector_variable_index( "ECM_attachment_point" );
	static int mechanical_strain_displacement_index = pCell->custom_data.find_vector_variable_index( "mechanical_strain_displacement" );
	pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value = pCell->custom_data.vector_variables[ECM_attachment_point_index].value;
	pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value -= pCell->position;
	pCell->custom_data[strain_index] = norm( pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value );

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

	if( strain_based_apoptosis( pCell ) )
	{
		pCell->phenotype.death.rates[apoptosis_index] = 9e9;
	}

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

	// plastoelastic mechanics
	static int spring_constant_index = pCell->custom_data.find_variable_index( "spring_constant" );
	static int relaxation_constant_index = pCell->custom_data.find_variable_index( "mechanical_relaxation_rate" );
	static int ECM_attachment_point_index = pCell->custom_data.find_vector_variable_index( "ECM_attachment_point" );
	static int mechanical_strain_displacement_index = pCell->custom_data.find_vector_variable_index( "mechanical_strain_displacement" );
	// first, update the cell's velocity based upon the elastic model
	axpy( &( pCell->velocity ) , pCell->custom_data[spring_constant_index] , pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value );

	// now, plastic mechanical relaxation
	static double plastic_temp_constant = -dt * pCell->custom_data[relaxation_constant_index];
	axpy( &(pCell->custom_data.vector_variables[ECM_attachment_point_index].value) , plastic_temp_constant , pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value );
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
	epithelium_submodel_info.microenvironment_variables.push_back( "TNF" );

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

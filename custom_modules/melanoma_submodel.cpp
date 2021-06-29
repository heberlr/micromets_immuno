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

double strain_based_proliferation( Cell* pCell )
{
	// std::vector<Cell*> neighbors = pCell->cells_in_my_container();//find cells in a neighbourhood of melanoma cells
	// static int lung_cell_type = get_cell_definition( "lung cell" ).type;
	// static int melanoma_cell_type = get_cell_definition( "melanoma cell" ).type;
	// int n = 0;
	// int CountNeighbors = 0;
	// Cell* pTestCell;
	// double radius_Cell = pCell->phenotype.geometry.radius;
	// while( n < neighbors.size() )
	// {
	// 	pTestCell = neighbors[n];
	// 	if ( pTestCell != pCell && pTestCell->phenotype.death.dead == false ){
	// 		double cell_cell_distance = sqrt((pTestCell->position[0]-pCell->position[0])*(pTestCell->position[0]-pCell->position[0])+(pTestCell->position[1]-pCell->position[1])*(pTestCell->position[1]-pCell->position[1]));
	// 		double radius_test_cell = pTestCell->phenotype.geometry.radius;
	// 		if( cell_cell_distance <= parameters.doubles("epsilon_distance")*(radius_Cell+radius_test_cell) ){
	// 			CountNeighbors++;
	// 		}
	// 	}
	// 	n++;
	// }
	//
	//
	// static double Xmin = microenvironment.mesh.bounding_box[0];
	// static double Ymin = microenvironment.mesh.bounding_box[1];
	// static double Xmax = microenvironment.mesh.bounding_box[3];
	// static double Ymax = microenvironment.mesh.bounding_box[4];
	//
	// if ( abs(pCell->position[0]-Xmin) < parameters.doubles("epsilon_distance")*radius_Cell ) return 0.0;
	// if ( abs(pCell->position[0]-Xmax) < parameters.doubles("epsilon_distance")*radius_Cell ) return 0.0;
	// if ( abs(pCell->position[1]-Ymin) < parameters.doubles("epsilon_distance")*radius_Cell ) return 0.0;
	// if ( abs(pCell->position[1]-Ymax) < parameters.doubles("epsilon_distance")*radius_Cell ) return 0.0;
	//
	// if ( CountNeighbors < 6 ) return 1.0;
	// else return 0.0;

	// static int strain_index = parameters.doubles("max_simple_pressure_TumorProl");//pCell->custom_data.find_variable_index( "mechanical_strain" );
	// static double  max_strain = parameters.doubles("max_mechanical_strain_TumorProl");
	//
	// if( pCell->custom_data[strain_index] < max_strain )
	// {
	// 	return pow( (max_strain - pCell->custom_data[strain_index])/max_strain, 1.0 );
	// }
	// return 0.0;

	static double  max_pressure = parameters.doubles("max_simple_pressure_TumorProl");
	if( pCell->state.simple_pressure < max_pressure )
	{
		return pow( (max_pressure - pCell->state.simple_pressure)/max_pressure, 1.0 );
	}
	return 0.0;
}

void melanoma_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	static int TNF_index = microenvironment.find_density_index( "TNF" );
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
	// update mechanical strain
	static int strain_index = pCell->custom_data.find_variable_index( "mechanical_strain" );
	static int ECM_attachment_point_index = pCell->custom_data.find_vector_variable_index( "ECM_attachment_point" );
	static int mechanical_strain_displacement_index = pCell->custom_data.find_vector_variable_index( "mechanical_strain_displacement" );
	pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value = pCell->custom_data.vector_variables[ECM_attachment_point_index].value;
	pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value -= pCell->position;
	pCell->custom_data[strain_index] = norm( pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value );

	//proliferation, mechanics, and chemokine secretion activated
	// Mechanical contribution to proliferation
	double mechanics_factor = strain_based_proliferation( pCell );
	//if (mechanics_factor == 0) pCell->phenotype.death.rates[apoptosis_index] = 9e9;

	// proliferation rate based on mechanical aspect
	int cycle_G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int cycle_S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
	pCell->phenotype.cycle.data.transition_rate(cycle_G0G1_index,cycle_S_index) = parameters.doubles("prolif_rate_CancerCell")*mechanics_factor;

	// pCell->custom_data["melanoma_cell_chemokine_secretion_activated"] = 1.0;
	// if( pCell->custom_data["melanoma_cell_chemokine_secretion_activated"] > 0.1 && phenotype.death.dead == false )
	// {
	// 	double rate = 1.0; //AV; // P;
	// 	rate /= pCell->custom_data["max_apoptosis_half_max"]; // Review this parameter, no make sense in melnoma dynamic
	// 	if( rate > 1.0 )
	// 	{ rate = 1.0; }
	// 	rate *= pCell->custom_data[ "melanoma_cell_chemokine_secretion_rate" ];
	//
	// 	phenotype.secretion.secretion_rates[chemokine_index] = rate;
	// 	phenotype.secretion.saturation_densities[chemokine_index] = 1.0;
	// }

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
	static int simple_pressure_index = pCell->custom_data.find_variable_index( "simple_pressure" );
	pCell->custom_data[simple_pressure_index]	= pCell->state.simple_pressure;
	//plastoelastic mechanics
	static int spring_constant_index = pCell->custom_data.find_variable_index( "spring_constant" );
	static int relaxation_constant_index = pCell->custom_data.find_variable_index( "mechanical_relaxation_rate" );
	static int ECM_attachment_point_index = pCell->custom_data.find_vector_variable_index( "ECM_attachment_point" );
	static int mechanical_strain_displacement_index = pCell->custom_data.find_vector_variable_index( "mechanical_strain_displacement" );

	pCell->custom_data[relaxation_constant_index] = 10e-3;

	// first, update the cell's velocity based upon the elastic model
	axpy( &( pCell->velocity ) , pCell->custom_data[spring_constant_index] , pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value );

	// now, plastic mechanical relaxation
	static double plastic_temp_constant = -dt * pCell->custom_data[relaxation_constant_index];
	axpy( &(pCell->custom_data.vector_variables[ECM_attachment_point_index].value) , plastic_temp_constant , pCell->custom_data.vector_variables[mechanical_strain_displacement_index].value );

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
	melanoma_submodel_info.microenvironment_variables.push_back( "TNF" );

	// what custom data do I need?
	//melanoma_submodel_info.cell_variables.push_back( "something" );
	// register the submodel
	melanoma_submodel_info.register_model();
	// set functions for the corresponding cell definition
	pCD = find_cell_definition( "melanoma cell" );
	pCD->functions.update_phenotype = melanoma_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = melanoma_submodel_info.mechanics_function;
	pCD->functions.contact_function = melanoma_contact_function;

	// Death rate for melanoma cells
	static int apoptosis_index = pCD->phenotype.death.find_death_model_index( "Apoptosis" );
	pCD->phenotype.death.rates[apoptosis_index] = parameters.doubles( "death_rate_CancerCell" );

	return;
}

void TCell_induced_apoptosis( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" );
	static int debris_index = microenvironment.find_density_index( "debris" );

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
		pCell->phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];

		pCell->functions.update_phenotype = NULL;
	}

	return;
}

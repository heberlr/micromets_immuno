#include "./immune_submodels.h"
#include <algorithm>

using namespace PhysiCell;

std::string immune_submodels_version = "0.4.0";
// Submodel_Information Immune_submodels_info; // not needed for now

Submodel_Information CD8_submodel_info;
Submodel_Information Macrophage_submodel_info;
Submodel_Information DC_submodel_info;
Submodel_Information CD4_submodel_info;

std::vector<Cell*> cells_to_move_from_edge;


std::vector<int> vascularized_voxel_indices;

// vector of valid positions
extern std::vector<std::vector <double>> valid_position;
// generator sample from distribution
extern std::default_random_engine generator;

extern AntigenLibrary AntigenLib;

// return true if out of bounds, within a tolerance
bool check_for_out_of_bounds( Cell* pC , double tolerance )
{
	static double Xmin = microenvironment.mesh.bounding_box[0];
	static double Ymin = microenvironment.mesh.bounding_box[1];
	static double Zmin = microenvironment.mesh.bounding_box[2];

	static double Xmax = microenvironment.mesh.bounding_box[3];
	static double Ymax = microenvironment.mesh.bounding_box[4];
	static double Zmax = microenvironment.mesh.bounding_box[5];

	static bool two_dimensions = default_microenvironment_options.simulate_2D;

	static bool setup_done = false;
	if( default_microenvironment_options.simulate_2D == true && setup_done == false )
	{
		Zmin = 0.0;
		Zmax = 0.0;
		setup_done = true;
	}

	if( pC->position[0] < Xmin + tolerance )
	{ return true; }
	if( pC->position[0] > Xmax - tolerance )
	{ return true; }

	if( pC->position[1] < Ymin + tolerance )
	{ return true; }
	if( pC->position[1] > Ymax - tolerance )
	{ return true; }

	if( two_dimensions )
	{ return false; }

	if( pC->position[2] < Zmin + tolerance )
	{ return true; }
	if( pC->position[2] > Zmax - tolerance )
	{ return true; }

	return false;
}

// return {push_x,push_y,push_z} of direction to nudge cell
std::vector<double> set_nudge_from_edge( Cell* pC , double tolerance )
{
	static double Xmin = microenvironment.mesh.bounding_box[0];
	static double Ymin = microenvironment.mesh.bounding_box[1];
	static double Zmin = microenvironment.mesh.bounding_box[2];

	static double Xmax = microenvironment.mesh.bounding_box[3];
	static double Ymax = microenvironment.mesh.bounding_box[4];
	static double Zmax = microenvironment.mesh.bounding_box[5];

	static bool two_dimensions = default_microenvironment_options.simulate_2D;

	static bool setup_done = false;
	if( default_microenvironment_options.simulate_2D == true && setup_done == false )
	{
		Zmin = 0.0;
		Zmax = 0.0;
		setup_done = true;
	}

	std::vector<double> nudge = {0,0,0};

	if( pC->position[0] < Xmin + tolerance )
	{ nudge[0] += 1; }
	if( pC->position[0] > Xmax - tolerance )
	{ nudge[0] -= 1; }

	if( pC->position[1] < Ymin + tolerance )
	{ nudge[1] += 1;  }
	if( pC->position[1] > Ymax - tolerance )
	{ nudge[1] -= 1;  }

	if( two_dimensions )
	{ normalize(nudge); return nudge; }

	if( pC->position[2] < Zmin + tolerance )
	{ nudge[2] += 1; }
	if( pC->position[2] > Zmax - tolerance )
	{ nudge[2] -= 1; }

	normalize(nudge);
	return nudge;
}

void nudge_out_of_bounds_cell( Cell* pC , double tolerance )
{
	std::vector<double> nudge = set_nudge_from_edge(pC,tolerance);

	// remove attachments
	pC->remove_all_attached_cells();

	// set velocity away rom edge
	pC->velocity = nudge;

	// set new position
	nudge *= tolerance;
	pC->position += nudge;

	// update in the data structure
	pC->update_voxel_in_container();

	// allow that cell to move and be movable
	pC->is_out_of_domain = false;
	pC->is_active = true;
	pC->is_movable= true;

	return;
}

void process_tagged_cells_on_edge( void )
{

	for( int n=0 ; n < cells_to_move_from_edge.size(); n++ )
	{
		Cell* pC = cells_to_move_from_edge[n];
		// std::cout << "moving cell from edge " << pC << " " << pC->type_name << std::endl;
		// replace_out_of_bounds_cell( cells_to_move_from_edge[n] , 10.0 );
		nudge_out_of_bounds_cell( pC , 10.0 );
	}
	//	if( cells_to_move_from_edge.size() > 0 )
	//	{ std::cout << std::endl; }

	return;
}

void choose_initialized_voxels( void )
{
	// read in percentage of tissue that's vascularised
	double percentage_vascularised = parameters.doubles("percentage_tissue_vascularized");
	int max_voxel_index = microenvironment.mesh.voxels.size() - 1;
	int number_of_vascularized_voxels = (int) ( percentage_vascularised/100.0 * ( max_voxel_index+1) );

	// choose which voxels are veins
	for( int n = 0 ; n < number_of_vascularized_voxels ; n++ )
	{
		int index_vascularised_voxel = (int) ( UniformRandom() * max_voxel_index );
		vascularized_voxel_indices.push_back( index_vascularised_voxel );
	}

	return;
}

void create_infiltrating_immune_cell( Cell_Definition* pCD )
{
	Cell* pC = create_cell( *pCD );

	std::vector<double> position = choose_vascularized_position();

	pC->assign_position( position );

	// Neoantigen signature to CD8 T cell - sample from antigen library (Heber)
	static int CD8_Tcell_type = get_cell_definition( "CD8 Tcell" ).type;
	static int neoantigen_signature_index = pC->custom_data.find_vector_variable_index( "neoantigen_signature" );
	if ( pC->type == CD8_Tcell_type ){
		pC->custom_data.vector_variables[neoantigen_signature_index].value = AntigenLib.sample_antigen();
	}

	return;
}

void create_infiltrating_immune_cell_initial( Cell_Definition* pCD )
{

	Cell* pC = create_cell( *pCD );
	// Sample from mesh positions (Heber)
	std::uniform_int_distribution<int> distribution_index(0, valid_position.size()-1);
	int index_sample = distribution_index(generator);
	pC->assign_position( valid_position[index_sample] );
	valid_position.erase (valid_position.begin()+index_sample);

	// randomly place cell intially
	/*double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Zmax = microenvironment.mesh.bounding_box[5];

	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0;
		Zmax = 0.0;
	}

	double Xrange = (Xmax - Xmin);
	double Yrange = (Ymax - Ymin);
	double Zrange = (Zmax - Zmin);

	// keep cells away from the outer edge

	Xmin += 0.1*Xrange;
	Ymin += 0.1*Yrange;
	Zmin = 0;

	Xrange *= 0.8;
	Yrange *= 0.8;
	Zrange = 0.0;

	// create some of each type of cell

	std::vector<double> position = {0,0,0};
	position[0] = Xmin + UniformRandom()*Xrange;
	position[1] = Ymin + UniformRandom()*Yrange;

	pC->assign_position( position );*/

	return;
}

std::vector<double> choose_vascularized_position( void )
{
	//extern std::vector<int> vascularized_voxel_indices;
	int my_voxel_index = (int) ( UniformRandom() * (vascularized_voxel_indices.size()-1) );
	int n = vascularized_voxel_indices[ my_voxel_index ] ;

	return microenvironment.mesh.voxels[n].center;
}

void create_infiltrating_immune_cell( std::string cell_name )
{
	create_infiltrating_immune_cell( find_cell_definition( cell_name ) );

	return;
}


void create_infiltrating_Tcell(void)
{
	static Cell_Definition* pCD = find_cell_definition( "CD8 Tcell" );
	create_infiltrating_immune_cell( pCD );

	return;
}

// (Adrianne) creating infiltrating CD4 cells
void create_infiltrating_CD4Tcell(void)
{
	static Cell_Definition* pCD = find_cell_definition( "CD4 Tcell" );
	create_infiltrating_immune_cell( pCD );

	return;
}

// (Adrianne) creating infiltrating DCs
void create_infiltrating_DC(void)
{
	static Cell_Definition* pCD = find_cell_definition( "DC" );
	create_infiltrating_immune_cell( pCD );

	return;
}

void create_infiltrating_macrophage(void)
{
	static Cell_Definition* pCD = find_cell_definition( "macrophage" );
	create_infiltrating_immune_cell( pCD );

	return;
}

void CD8_Tcell_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2 , double dt )
{
	//std::cout << pC1 << " " << pC1->type_name << " contact with " << pC2 << " " << pC2->type_name << std::endl;
	// elastic adhesions
	standard_elastic_contact_function( pC1,p1, pC2, p2, dt );

	// increase contact time of cell you are attacking
	#pragma omp critical
	{ pC2->custom_data["TCell_contact_time"] += dt; /*std::cout << "TCell_contact_time: " << pC2->custom_data["TCell_contact_time"] << std::endl;*/}

	return;
}

void CD8_Tcell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// MISSING CHANGE DEATH RATE (DEPEND ON GENERATION)

	static int debris_index = microenvironment.find_density_index( "debris");

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		pCell->functions.custom_cell_rule = NULL;

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
		return;
	}

	return;
}

void CD8_Tcell_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		pCell->functions.custom_cell_rule = NULL;

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
		return;
	}

	// bounds check
	if( check_for_out_of_bounds( pCell , 10.0 ) )
	{
		#pragma omp critical
		{ cells_to_move_from_edge.push_back( pCell ); }
		// replace_out_of_bounds_cell( pCell, 10.0 );
		// return;
	}

	//if (pCell->state.neighbors.size()>0){ std::cout << "ID: " << pCell->ID << " Neighbors: " << pCell->state.neighbors.size() << std::endl; exit(0);}

	// if I am not adhered to a cell, turn motility on
	if( pCell->state.neighbors.size() == 0 )
	{ phenotype.motility.is_motile = true; }
	else
	{ phenotype.motility.is_motile = false; }

	// check for contact with melanoma cell

	// if I'm adhered to something ...
	if( pCell->state.number_of_attached_cells() > 0 ) // pCell->state.neighbors.size() > 0 )
	{
		// decide whether to detach
		bool detach_me = false;

		if( UniformRandom() < dt / ( pCell->custom_data["cell_attachment_lifetime"] + 1e-15 ) )
		{ detach_me = true; }

		// if I detach, go through the process
		if( detach_me )
		{
			pCell->remove_all_attached_cells();
			// resume motile behavior
			phenotype.motility.is_motile = true;
		}
		return;
	}

	// I'm not attached, look for cells nearby and try to attach

	// if this returns non-NULL, we're now attached to a cell
	if( immune_cell_check_neighbors_for_attachment( pCell , dt) )
	{
		// set motility off
		phenotype.motility.is_motile = false;
		return;
	}

	return;
}

void immune_cell_motility_direction( Cell* pCell, Phenotype& phenotype , double dt )
{
	if( phenotype.death.dead == true )
	{
		phenotype.motility.migration_speed = 0.0;
		return;
	}

	static int TNF_index = microenvironment.find_density_index( "TNF");
	static int debris_index = microenvironment.find_density_index( "debris");

	// if not activated, chemotaxis along debris
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(debris_index);
	normalize( &phenotype.motility.migration_bias_direction );
	if( pCell->custom_data["activated_immune_cell"] < 0.5 )
	{ return; }

	// if activated, follow the weighted direction

	phenotype.motility.migration_bias_direction *= pCell->custom_data["sensitivity_to_debris_chemotaxis"];

	std::vector<double> gradC = pCell->nearest_gradient(TNF_index);
	normalize( &gradC );
	gradC *= pCell->custom_data["sensitivity_to_TNF_chemotaxis"];

	phenotype.motility.migration_bias_direction += gradC;

	normalize( &( phenotype.motility.migration_bias_direction) );


	return;
}

void macrophage_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int apoptosis_index = phenotype.death.find_death_model_index( "Apoptosis" );
	static Cell_Definition* pCD = find_cell_definition( "macrophage" );
	static int TNF_index = microenvironment.find_density_index( "TNF");
	static int debris_index = microenvironment.find_density_index( "debris");

	// no apoptosis until activation (resident macrophages in constant number for homeostasis)
	if( pCell->custom_data["activated_immune_cell"] < 0.5 )
	{ phenotype.death.rates[apoptosis_index] = 0.0; }
	else
	{ phenotype.death.rates[apoptosis_index] = pCD->phenotype.death.rates[apoptosis_index]; }

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		pCell->functions.custom_cell_rule = NULL;

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
		return;
	}

	// check for cells to eat
	std::vector<Cell*> neighbors = pCell->cells_in_my_container();

	// at least one of the cells is pCell
	if( neighbors.size() < 2 )
	{ return; }

	// (Adrianne) get type of CD8+ T cell and CD4+ t CELL
	static int CD8_Tcell_type = get_cell_definition( "CD8 Tcell" ).type;
	static int CD4_Tcell_type = get_cell_definition( "CD4 Tcell" ).type;

	// (Adrianne) if there is a T cell in a mac's neighbourhood AND a mac has already begin phagocytosing, then there will be changes to the macs actions
	int n = 0;
	Cell* pContactCell = neighbors[n];
	while( n < neighbors.size() )
	{
		pContactCell = neighbors[n];

		double cell_cell_distance = sqrt((pContactCell->position[0]-pCell->position[0])*(pContactCell->position[0]-pCell->position[0])+(pContactCell->position[1]-pCell->position[1])*(pContactCell->position[1]-pCell->position[1]));
		double radius_mac = pCell->phenotype.geometry.radius; // (Adrianne) radius of DC)
		double radius_test_cell = pContactCell->phenotype.geometry.radius; // (Adrianne) radius of test cell)

		// (Adrianne) if it is not me, not dead and is a CD8 T cell that is within a very short distance from me, I will stop secreting TNF
		if( pContactCell != pCell && pContactCell->phenotype.death.dead == false && pContactCell->type == CD8_Tcell_type
			&& pCell->custom_data["activated_immune_cell"] > 0.5 && cell_cell_distance<=parameters.doubles("epsilon_distance")*(radius_mac+radius_test_cell))
			{
				phenotype.secretion.secretion_rates[TNF_index] = 0;// (Adrianne) contact with CD8 T cell turns off TNF secretion
				n=neighbors.size();
			}
			// (Adrianne) if it is not me, not dead and is a CD4 T cell that is within a very short distance from me, I will be able to phagocytose tumor (but not neccesarily dead) cells
		else if( pContactCell != pCell && pContactCell->phenotype.death.dead == false && pContactCell->type == CD4_Tcell_type
			&& pCell->custom_data["activated_immune_cell"] > 0.5 && cell_cell_distance<=parameters.doubles("epsilon_distance")*(radius_mac+radius_test_cell))
			{
				pCell->custom_data["ability_to_phagocytose_melanoma_cell"] = 1; // (Adrianne) contact with CD4 T cell induces macrophage's ability to phagocytose tumor cells
				n=neighbors.size();
			}
			n++;
	}

	// (Adrianne) if macrophage volume exceeds a threshold value we say it is "exhausted" and unable to phagocytose until it's volume drops below this threshold
	if( pCell->phenotype.volume.total> pCell->custom_data["threshold_macrophage_volume"])
	{
		// (Adrianne) when a macrophage is in an exhausted state it has a death rate  2.1e-4
		phenotype.death.rates[apoptosis_index] = pCell->custom_data["exhausted_macrophage_death_rate"];
		return;
	}

	// (Adrianne) obtain index for tracking time when next phagocytosis event is possible
	int time_to_next_phagocytosis_index = pCell->custom_data.find_variable_index( "time_to_next_phagocytosis" );
	// (Adrianne) check if still phagocytosing something, added if statement to say that if cell is still internalising current material not to phagocytose anything else
	if( pCell->custom_data.variables[time_to_next_phagocytosis_index].value>PhysiCell_globals.current_time )
	{return;}

	double probability_of_phagocytosis = pCell->custom_data["phagocytosis_rate"] * dt;

	// (Adrianne) add an additional variable that is the time taken to ingest material
	double material_internalisation_rate = pCell->custom_data["material_internalisation_rate"];

	n = 0;
	Cell* pTestCell = neighbors[n];
	while( n < neighbors.size() )
	{
		pTestCell = neighbors[n];
		// if it is not me and not a macrophage
		if( pTestCell != pCell && pTestCell->phenotype.death.dead == true &&
			UniformRandom() < probability_of_phagocytosis ) // && // remove in v 3.2
			//			pTestCell->phenotype.volume.total < max_phagocytosis_volume ) / remove in v 3.2
			{
				{
					// (Adrianne) obtain volume of cell to be ingested
					double volume_ingested_cell = pTestCell->phenotype.volume.total;

					pCell->ingest_cell( pTestCell );

					// (Adrianne) macrophage cannot phagocytose again until it has elapsed the time taken to phagocytose the material
					double time_to_ingest = volume_ingested_cell*material_internalisation_rate;// convert volume to time taken to phagocytose
					// (Adrianne) update internal time vector in macrophages that tracks time it will spend phagocytosing the material so they can't phagocytose again until this time has elapsed
					pCell->custom_data.variables[time_to_next_phagocytosis_index].value = PhysiCell_globals.current_time+time_to_ingest;
				}

				// activate the cell
				phenotype.secretion.secretion_rates[TNF_index] =
				pCell->custom_data["activated_TNF_secretion_rate"]; // 10;
				phenotype.secretion.saturation_densities[TNF_index] = 1;

				phenotype.secretion.uptake_rates[TNF_index] = 0.0;

				phenotype.motility.migration_speed = pCell->custom_data["activated_speed"];

				//adding TNF uptake by phagocytes (Heber)

				pCell->custom_data["activated_immune_cell"] = 1.0;

				return;
			}
			else if( pTestCell != pCell && pCell->custom_data["ability_to_phagocytose_melanoma_cell"]== 1 && /*pTestCell->custom_data[nP]>1 &&*/
			UniformRandom() < probability_of_phagocytosis ) // (Adrianne) macrophages that have been activated by T cells can phagocytose infected cells that contain at least 1 pathogen protein
			{
				{
					// (Adrianne) obtain volume of cell to be ingested
					double volume_ingested_cell = pTestCell->phenotype.volume.total;

					pCell->ingest_cell( pTestCell );

					// (Adrianne) macrophage cannot phagocytose again until it has elapsed the time taken to phagocytose the material
					double time_to_ingest = volume_ingested_cell*material_internalisation_rate;// convert volume to time taken to phagocytose
					// (Adrianne) update internal time vector in macrophages that tracks time it will spend phagocytosing the material so they can't phagocytose again until this time has elapsed
					pCell->custom_data.variables[time_to_next_phagocytosis_index].value = PhysiCell_globals.current_time+time_to_ingest;
				}

				// activate the cell
				phenotype.secretion.secretion_rates[TNF_index] =
				pCell->custom_data["activated_TNF_secretion_rate"]; // 10;
				phenotype.secretion.saturation_densities[TNF_index] = 1;

				phenotype.secretion.uptake_rates[TNF_index] = 0.0;

				phenotype.motility.migration_speed = pCell->custom_data["activated_speed"];

				pCell->custom_data["activated_immune_cell"] = 1.0;

				return;
			}
			n++;
		}
	return;
}

void macrophage_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
				static int debris_index = microenvironment.find_density_index( "debris");

				if( phenotype.death.dead == true )
				{
					pCell->functions.update_phenotype = NULL;
					pCell->functions.custom_cell_rule = NULL;

					phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
					return;
				}

				// bounds check
				if( check_for_out_of_bounds( pCell , 10.0 ) )
				{
					#pragma omp critical
					{ cells_to_move_from_edge.push_back( pCell ); }
					// replace_out_of_bounds_cell( pCell, 10.0 );
					// return;
				}

				//	// death check
				//	if( phenotype.death.dead == true )
				//	{ remove_all_adhesions( pCell ); }

				return;
}

void DC_contact_function( Cell* pC1, Phenotype& p1, Cell* pC2, Phenotype& p2 , double dt )
{
	//std::cout << pC1 << " " << pC1->type_name << " contact with " << pC2 << " " << pC2->type_name << std::endl;
	// elastic adhesions
	standard_elastic_contact_function( pC1,p1, pC2, p2, dt );

	return;
}
// (Adrianne) DC phenotype function
void DC_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	// (Adrianne) get type of CD8+ T cell
	static int CD8_Tcell_type = get_cell_definition( "CD8 Tcell" ).type;
	Cell* pTempCell = NULL;

	// Index neoantigen signature
	static int neoantigen_signature_index = pCell->custom_data.find_vector_variable_index( "neoantigen_signature" );

	// (Adrianne) if DC is already activated, then check whether it leaves the tissue
	if( pCell->custom_data["activated_immune_cell"] >  0.5 && UniformRandom() < 0.002)
	{
		extern double DM; //declare existance of DC lymph
		// (Adrianne) DC leaves the tissue and so we lyse that DC
		std::cout<<"DC leaves tissue"<<std::endl;
		pCell->lyse_cell();
		#pragma omp critical
		{
			DM++; // add one
		}
		return;

	}
	else if( pCell->custom_data["activated_immune_cell"] > 0.5 ) // (Adrianne) activated DCs that don't leave the tissue can further activate CD8s increasing their proliferation rate and attachment rates
	{

		std::vector<Cell*> neighbors = pCell->cells_in_my_container(); // (Adrianne) find cells in a neighbourhood of DCs
		int n = 0;
		Cell* pTestCell = neighbors[n];
		while( n < neighbors.size() )
		{
			if ( pTestCell != pCell && pTestCell->phenotype.death.dead == false && pTestCell->type == CD8_Tcell_type)
			pTestCell = neighbors[n];

			// (Adrianne) find the euclidean distance between the DC and the cell it's testing
			double cell_cell_distance = sqrt((pTestCell->position[0]-pCell->position[0])*(pTestCell->position[0]-pCell->position[0])+(pTestCell->position[1]-pCell->position[1])*(pTestCell->position[1]-pCell->position[1]));
			double radius_DC = pCell->phenotype.geometry.radius; // (Adrianne) radius of DC)
			double radius_test_cell = pTestCell->phenotype.geometry.radius; // (Adrianne) radius of test cell)

			// (Adrianne) check if any neighbour cells are live T cells and that they are close enough to the DC
			if( pTestCell != pCell && pTestCell->phenotype.death.dead == false && pTestCell->type == CD8_Tcell_type
				&& cell_cell_distance<=parameters.doubles("epsilon_distance")*(radius_DC+radius_test_cell))
				{

					pTestCell-> custom_data["cell_attachment_rate"] = parameters.doubles("DC_induced_CD8_attachment"); // (Adrianne) DC induced T cell attachement rate

					// (Adrianne) finding the G0G1 and S phase index and setting the transition rate to be non zero so that CD8 T cells start proliferating after interacting with DC
					int cycle_G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
					int cycle_S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
					pTestCell->phenotype.cycle.data.transition_rate(cycle_G0G1_index,cycle_S_index) = parameters.doubles("DC_induced_CD8_proliferation");

					n = neighbors.size();

				}

				n++;

			}

			return;
		}
		else
		{
			//Activation of DC cells - attach with death melanoma cells
			// if this returns non-NULL, we're now attached to at least one cell
			if( immune_cell_check_neighbors_for_attachment( pCell , dt) )
			{
				// Look around by melanoma cells attached)
				for (int i = 0; i < pCell->state.number_of_attached_cells(); i++){
					pTempCell = pCell->state.attached_cells[i];
					//Add antigen to library
					#pragma omp critical
					{
						AntigenLib.add_antigen(pTempCell->custom_data.vector_variables[neoantigen_signature_index].value, PhysiCell_globals.current_time );
						//AntigenLib.print();
						pCell->custom_data["activated_immune_cell"] = 1.0;
						phenotype.motility.is_motile = false;
					}
				}
			}
		}
		return;
	}

// (Adrianne) DC mechanics function
void DC_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		pCell->functions.custom_cell_rule = NULL;

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
		return;
	}

	// bounds check
	if( check_for_out_of_bounds( pCell , 10.0 ) )
	{
		#pragma omp critical
		{ cells_to_move_from_edge.push_back( pCell ); }
		// replace_out_of_bounds_cell( pCell, 10.0 );
		// return;
	}

	//	// death check
	//	if( phenotype.death.dead == true )
	//	{ remove_all_adhesions( pCell ); }

	return;
}

// (Adrianne CD4 phenotype function
void CD4_Tcell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
	int cycle_G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int cycle_S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );

	static int apoptosis_index = pCell->phenotype.death.find_death_model_index( "apoptosis" );

	if(pCell->phenotype.cycle.data.elapsed_time_in_phase<6 &&  pCell->phenotype.cycle.data.current_phase_index==0)
	{
		pCell->phenotype.death.rates[apoptosis_index] = 100; // new death rate of T cells when they have exceeded generation

		pCell->phenotype.cycle.data.transition_rate(cycle_G0G1_index,cycle_S_index) = 0;

	}
	return;
}

// (Adrianne) CD4 mechanics function
void CD4_Tcell_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");

	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		pCell->functions.custom_cell_rule = NULL;

		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"];
		return;
	}
	// bounds check
	if( check_for_out_of_bounds( pCell , 10.0 ) )
	{
		#pragma omp critical
		{ cells_to_move_from_edge.push_back( pCell ); }
		// replace_out_of_bounds_cell( pCell, 10.0 );
		// return;
	}

//	// death check
//	if( phenotype.death.dead == true )
//	{ remove_all_adhesions( pCell ); }

	return;
}

void immune_submodels_setup( void )
{
	Cell_Definition* pCD;

	//
	// set up CD8 Tcells
	// set version info
	CD8_submodel_info.name = "CD8 Tcell model";
	CD8_submodel_info.version = immune_submodels_version;
	// set functions
	CD8_submodel_info.main_function = NULL;
	CD8_submodel_info.phenotype_function = CD8_Tcell_phenotype;
	CD8_submodel_info.mechanics_function = CD8_Tcell_mechanics;
	// what microenvironment variables do you expect?
	CD8_submodel_info.microenvironment_variables.push_back( "TNF" );
	// what custom data do I need?
	//CD8_submodel_info.cell_variables.push_back( "something" );
	// register the submodel
	CD8_submodel_info.register_model();
	// set functions for the corresponding cell definition
	pCD = find_cell_definition( "CD8 Tcell" );
	pCD->functions.update_phenotype = CD8_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = CD8_submodel_info.mechanics_function;
	pCD->functions.contact_function = CD8_Tcell_contact_function;

	// set up macrophages
	Macrophage_submodel_info = CD8_submodel_info; // much shared information
	// set version info
	Macrophage_submodel_info.name = "macrophage model";
	Macrophage_submodel_info.version = immune_submodels_version;
	// set functions
	Macrophage_submodel_info.main_function = NULL;
	Macrophage_submodel_info.phenotype_function = macrophage_phenotype;
	Macrophage_submodel_info.mechanics_function = macrophage_mechanics;
	// what microenvironment variables do you expect?
	// nothing unique
	// what custom data do I need?
	//CD8_submodel_info.cell_variables.push_back( "something" );
	// register the submodel
	Macrophage_submodel_info.register_model();
	// set functions for the corresponding cell definition
	pCD = find_cell_definition( "macrophage" );
	pCD->functions.update_phenotype = Macrophage_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = Macrophage_submodel_info.mechanics_function;
	pCD->functions.update_migration_bias = immune_cell_motility_direction;


	// (Adrianne) set up DC submodel info
	DC_submodel_info = CD8_submodel_info; // much shared information
	// set version info
	DC_submodel_info.name = "DC model";
	DC_submodel_info.version = immune_submodels_version;
	// set functions
	DC_submodel_info.main_function = NULL;
	DC_submodel_info.phenotype_function = DC_phenotype;
	DC_submodel_info.mechanics_function = DC_mechanics;
	// what microenvironment variables do you expect?
	// nothing unique
	// what custom data do I need?
	//CD8_submodel_info.cell_variables.push_back( "something" );
	// register the submodel
	DC_submodel_info.register_model();
	// set functions for the corresponding cell definition
	pCD = find_cell_definition( "DC" );
	pCD->functions.update_phenotype = DC_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = DC_submodel_info.mechanics_function;
	pCD->functions.contact_function = DC_contact_function; // (Heber) Attached cells
	pCD->functions.update_migration_bias = immune_cell_motility_direction;

	// (Adrianne) set up CD4 Tcells ** we don't want CD4's to do anything expect be recruited to the tissue and migrate in tissue
	// set version info
	CD4_submodel_info = CD8_submodel_info; // much shared information
	CD4_submodel_info.name = "CD4 Tcell model";
	CD4_submodel_info.version = immune_submodels_version;

	// set functions
	CD4_submodel_info.main_function = NULL;
	CD4_submodel_info.phenotype_function = CD4_Tcell_phenotype;
	CD4_submodel_info.mechanics_function = CD4_Tcell_mechanics;
	// what microenvironment variables do you expect?
	CD4_submodel_info.microenvironment_variables.push_back( "TNF" );
	// what custom data do I need?
	//CD8_submodel_info.cell_variables.push_back( "something" );
	// register the submodel
	CD4_submodel_info.register_model();
	// set functions for the corresponding cell definition
	pCD = find_cell_definition( "CD4 Tcell" );
	pCD->functions.update_phenotype = CD4_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = CD4_submodel_info.mechanics_function;
}

Cell* check_for_live_neighbor_for_interaction( Cell* pAttacker , double dt )
{
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container();
	int i = 0;
	while( i < nearby.size() )
	{
		// don't try to kill yourself
		if( nearby[i] != pAttacker && nearby[i]->phenotype.death.dead == false )
		{ return nearby[i]; }
		i++;
	}
	return NULL;
}

Cell* check_for_dead_neighbor_for_interaction( Cell* pAttacker , double dt )
{
std::vector<Cell*> nearby = pAttacker->cells_in_my_container();
int i = 0;
while( i < nearby.size() )
{
	// don't try to kill yourself
	if( nearby[i] != pAttacker && nearby[i]->phenotype.death.dead == true )
	{ return nearby[i]; }
	i++;
}
return NULL;
}

bool attempt_immune_cell_attachment( Cell* pAttacker, Cell* pTarget , double dt )
{
static int CD8_Tcell_type = get_cell_definition( "CD8 Tcell" ).type;
static int DC_type = get_cell_definition( "DC" ).type;
static int neoantigen_signature_index = pAttacker->custom_data.find_vector_variable_index( "neoantigen_signature" );

// if the target is not melanoma cell, give up for CD8 attack
static int melanoma_type = get_cell_definition( "melanoma cell" ).type;
if ( pTarget->type != melanoma_type && pAttacker->type == CD8_Tcell_type)
{ return false; }

// if the target is not melanoma or lung cell, give up for DC attack
static int lung_type = get_cell_definition( "lung cell" ).type;
if ( pTarget->type != melanoma_type && pTarget->type != lung_type  && pAttacker->type == DC_type)
{ return false; }

// if the target cell is dead, give up for CD8 attack
if( pTarget->phenotype.death.dead == true && pAttacker->type == CD8_Tcell_type)
{ return false; }

// if the target cell is not dead, give up for DC attach (melanoma dead cell presenting antigen)
if( pTarget->phenotype.death.dead != true && pAttacker->type == DC_type )
{ return false; }

// if the target cell is too far away, give up
std::vector<double> displacement = pTarget->position - pAttacker->position;
double distance_scale = norm( displacement );

// better: use mechanics constants
if( distance_scale > pAttacker->custom_data["max_attachment_distance"] )
{ return false; }

// now, get the attachment probability

double attachment_probability = pAttacker->custom_data["cell_attachment_rate"] * dt;

// neoantigen recognition for CD8 T cell attacking
if ( pAttacker->type == CD8_Tcell_type){
	int HammingDistance = Hamming_Distance(pAttacker->custom_data.vector_variables[neoantigen_signature_index].value, pTarget->custom_data.vector_variables[neoantigen_signature_index].value);
	attachment_probability *= (1.0/(1.0 + pow(HammingDistance/parameters.doubles( "half_hamming_distance" ), parameters.doubles( "hill_coefficient_recognition" )) ) );
}
// don't need to cap it at 1.00: if prob > 100%,
// then this statement always evaluates as true,
// just the same as capping probability at 100%
if( UniformRandom() <= attachment_probability )
{
	attach_cells( pAttacker, pTarget );
	return true;
}

return false;
}

Cell* immune_cell_check_neighbors_for_attachment( Cell* pAttacker , double dt )
{
	static int CD8_Tcell_type = get_cell_definition( "CD8 Tcell" ).type;
	std::vector<Cell*> nearby = pAttacker->cells_in_my_container();
	int i = 0;
	while( i < nearby.size() )
	{
		// don't try to kill yourself
		if( nearby[i] != pAttacker )
		{
			// CD8 T cell attaches only to one cell
			if( attempt_immune_cell_attachment( pAttacker, nearby[i] , dt ) && pAttacker->type == CD8_Tcell_type )
			{
				return nearby[i];
			}
		}
		i++;
	}

	// Cell has more than one attached cell (Heber)
	if ( pAttacker->state.number_of_attached_cells() > 0 ) return pAttacker->state.attached_cells[0];

	return NULL;
}

int recruited_Tcells = 0;
int recruited_macrophages = 0;
int recruited_CD4Tcells = 0;
int recruited_DCs = 0;

double first_macrophage_recruitment_time = 9e9;
double first_CD8_T_cell_recruitment_time = 9e9;
double first_CD4_T_cell_recruitment_time = 9e9;
double first_DC_recruitment_time = 9e9;

void immune_cell_recruitment( double dt )
{
	static int TNF_index =
	microenvironment.find_density_index("TNF");

	static double dt_immune = parameters.doubles( "immune_dt" );
	static double t_immune = 0.0;
	static double t_last_immune = 0.0;
	static double t_next_immune = 0.0;

	static double tolerance = 0.1 * diffusion_dt;

	// is it time for the next immune recruitment?
	if( t_immune > t_next_immune- tolerance )
	{
		double elapsed_time = (t_immune - t_last_immune );

		// macrophage recruitment

		static double macrophage_recruitment_rate = parameters.doubles( "macrophage_max_recruitment_rate" );
		static double M_min_signal = parameters.doubles( "macrophage_recruitment_min_signal" );
		static double M_sat_signal = parameters.doubles( "macrophage_recruitment_saturation_signal" );
		static double M_max_minus_min = M_sat_signal - M_min_signal;

		double total_rate = 0;
		// integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV
		double total_scaled_signal= 0.0;
		for( int n=0; n<microenvironment.mesh.voxels.size(); n++ )
		{
			// (signal(x)-signal_min)/(signal_max/signal_min)
			double dRate = ( microenvironment(n)[TNF_index] - M_min_signal );
			dRate /= M_max_minus_min;
			// crop to [0,1]
			if( dRate > 1 )
			{ dRate = 1; }
			if( dRate < 0 )
			{ dRate = 0; }
			total_rate += dRate;
		}
		// multiply by dV and rate_max
		total_scaled_signal = total_rate;

		total_rate *= microenvironment.mesh.dV;
		total_rate *= macrophage_recruitment_rate;

		// expected number of new macrophages
		double number_of_new_cells_prob = total_rate * elapsed_time ;

		int number_of_new_cells_int = floor( number_of_new_cells_prob );
		double alpha = number_of_new_cells_prob - number_of_new_cells_int;

		//STOCHASTIC PORTION

		if( UniformRandom()< alpha )
		{
			number_of_new_cells_int++;
		}
		recruited_macrophages += number_of_new_cells_int;

		if( number_of_new_cells_int )
		{
			if( t_immune < first_macrophage_recruitment_time )
			{ first_macrophage_recruitment_time = t_immune; }

			std::cout << "\tRecruiting " << number_of_new_cells_int << " macrophages ... " << std::endl;

			for( int n = 0; n < number_of_new_cells_int ; n++ )
			{ create_infiltrating_macrophage(); }
		}

		// CD8 Tcell recruitment (Michael) changed to take floor of ODE value
		extern double TCt;
		extern std::vector<int>historyTc;

		int number_of_new_cells = (int) floor( TCt );
		TCt -= number_of_new_cells;

		std::rotate(historyTc.rbegin(),historyTc.rbegin()+1,historyTc.rend());
		historyTc.front() = number_of_new_cells;

		recruited_Tcells += historyTc.back();

		if( historyTc.back() )
		{
			if( t_immune < first_CD8_T_cell_recruitment_time )
			{ first_CD8_T_cell_recruitment_time = t_immune; }

			std::cout << "\tRecruiting " << historyTc.back() << " CD8 T cells ... " << std::endl;

			for( int n = 0; n < historyTc.back(); n++ )
			{ create_infiltrating_Tcell(); }
		}

		// CD4 recruitment (Michael) changed to take floor of ODE value
		extern double Tht;
		extern std::vector<int>historyTh;

		number_of_new_cells = (int) floor( Tht );
		Tht-=number_of_new_cells;

		std::rotate(historyTh.rbegin(),historyTh.rbegin()+1,historyTh.rend());
		historyTh.front() = number_of_new_cells;

		recruited_CD4Tcells += historyTh.back();

		if( historyTh.back() )
		{
			if( t_immune < first_CD4_T_cell_recruitment_time )
			{ first_CD4_T_cell_recruitment_time = t_immune; }

			std::cout << "\tRecruiting " << historyTh.back() << " CD4 T cells ... " << std::endl;

			for( int n = 0; n < historyTh.back() ; n++ )
			{ create_infiltrating_CD4Tcell(); }
		}

		// (Adrianne) DC recruitment - *** This section will be changed to be Tarun's model  so I've left recruitment parameters to be mac cell parameters**
		static double DC_recruitment_rate = parameters.doubles( "DC_max_recruitment_rate" );
		static double DC_min_signal = parameters.doubles( "DC_recruitment_min_signal" );
		static double DC_sat_signal = parameters.doubles( "DC_recruitment_saturation_signal" );
		static double DC_max_minus_min = DC_sat_signal - DC_min_signal;

		total_rate = 0;
		// integrate \int_domain r_max * (signal-signal_min)/(signal_max-signal_min) * dV
		total_scaled_signal= 0.0;
		for( int n=0; n<microenvironment.mesh.voxels.size(); n++ )
		{
			// (signal(x)-signal_min)/(signal_max/signal_min)
			double dRate = ( microenvironment(n)[TNF_index] - DC_min_signal );
			dRate /= DC_max_minus_min;
			// crop to [0,1]
			if( dRate > 1 )
			{ dRate = 1; }
			if( dRate < 0 )
			{ dRate = 0; }
			total_rate += dRate;
		}
		// multiply by dV and rate_max
		total_scaled_signal = total_rate;

		total_rate *= microenvironment.mesh.dV;
		total_rate *= DC_recruitment_rate;

		// expected number of new DCs
		number_of_new_cells_prob = total_rate * elapsed_time ;

		number_of_new_cells_int = floor( number_of_new_cells_prob );
		alpha = number_of_new_cells_prob - number_of_new_cells_int;

		//STOCHASTIC PORTION

		if( UniformRandom()< alpha )
		{
			number_of_new_cells_int++;
		}
		recruited_DCs += number_of_new_cells_int;

		if( number_of_new_cells_int )
		{
			if( t_immune < first_DC_recruitment_time )
			{ first_DC_recruitment_time = t_immune; }

			std::cout << "\tRecruiting " << number_of_new_cells_int << " DCs ... " << std::endl;

			for( int n = 0; n < number_of_new_cells_int ; n++ )
			{ create_infiltrating_DC(); }
		}


		t_last_immune = t_immune;
		t_next_immune = t_immune + dt_immune;

	}


	t_immune += dt;

	return;
}

void initial_immune_cell_placement( void )
{
	Cell_Definition* pCD8 = find_cell_definition( "CD8 Tcell" );
	Cell_Definition* pMF = find_cell_definition( "macrophage" );
	Cell_Definition* pDC = find_cell_definition( "DC" );
	Cell_Definition* pCD4 = find_cell_definition( "CD4 Tcell" );

	// CD8+ T cells;
	for( int n = 0 ; n < parameters.ints("number_of_CD8_Tcells") ; n++ )
	{ create_infiltrating_immune_cell( pCD8 ); }

	// macrophages
	for( int n = 0 ; n < parameters.ints("number_of_macrophages") ; n++ )
	{ create_infiltrating_immune_cell_initial( pMF ); }

	// DC
	for( int n = 0 ; n < parameters.ints("number_of_DCs") ; n++ )
	{ create_infiltrating_immune_cell_initial( pDC ); }

	// CD4+ T cells
	for( int n = 0 ; n < parameters.ints("number_of_CD4_Tcells") ; n++ )
	{ create_infiltrating_immune_cell_initial( pCD4 ); }

	return;
}

int Hamming_Distance(const std::vector<double> &Seq1, const std::vector<double> &Seq2){
  int count = 0;
  for (int i=0; i < Seq1.size(); i++){
    if (Seq1[i] != Seq2[i]) count++;
  }
  return count;
}

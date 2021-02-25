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
    // Heber - Concentration of danger signals on cell
    static int danger_signals_index = microenvironment.find_density_index( "danger signals" );
	
	// receptor dynamics 
	// requires faster time scale - done in main function 
	
	// pathogen dynamics model 
	//internal_pathogen_dynamics_info.phenotype_function(pCell,phenotype,dt); 
	// internal_pathogen_model(pCell,phenotype,dt);
	
	// pathogen response model 
	//internal_pathogen_response_model_info.phenotype_function(pCell,phenotype,dt); 
	// internal_pathogen_response_model(pCell,phenotype,dt);

    //After is possible integrate a intracellular model to danger signals
	pCell->custom_data["danger_signals_intracellular"] =  pCell->phenotype.secretion.secretion_rates[danger_signals_index];//pCell->phenotype.molecular.internalized_total_substrates[ danger_signals_index ];
    phenotype.secretion.saturation_densities[danger_signals_index] = 1.0;
	// T-cell based death
	TCell_induced_apoptosis(pCell, phenotype, dt ); 
	
	// if I am dead, remove all adhesions 
	if( phenotype.death.dead == true )
	{
		// detach all attached cells 
		// remove_all_adhesions( pCell ); 
		
		phenotype.secretion.secretion_rates[debris_index] = pCell->custom_data["debris_secretion_rate"]; 
	}
	
	// cell secretion belongs in pathogen response 
	
	// if I am dead, make sure to still secrete the chemokine 
	static int chemokine_index = microenvironment.find_density_index( "chemokine" ); 
    
    
    
    // Heber - Cell mutation
    /*double probability_of_mutation = pCell->custom_data["mutation_rate"] * dt;
    if (pCell->phenotype.molecular.internalized_total_substrates[ danger_signals_index ] == 0 && UniformRandom() < probability_of_mutation)  
        phenotype.secretion.secretion_rates[danger_signals_index] = 1.0;*/
            
	// Heber - Release chemokine if danger signals is more than a certain threshold 
    //if( (pCell->nearest_density_vector())[ danger_signals_index ] > 0.5 )  std::cout << "TEST: " << (pCell->nearest_density_vector())[ danger_signals_index ] << std::endl;
    if( (pCell->nearest_density_vector())[ danger_signals_index ] >= 0.8 ) 
	{
		pCell->custom_data["mutated_cell_chemokine_secretion_activated"] = 1.0; 
	}
	if( pCell->custom_data["mutated_cell_chemokine_secretion_activated"] > 0.1 && phenotype.death.dead == false )
	{
		double rate = 1.0; //AV; // P; 
		rate /= pCell->custom_data["max_apoptosis_half_max"];
		if( rate > 1.0 )
		{ rate = 1.0; }
		rate *= pCell->custom_data[ "mutated_cell_chemokine_secretion_rate" ];

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

void epithelium_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	static int debris_index = microenvironment.find_density_index( "debris");
	
	pCell->is_movable = false; 
	
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
	epithelium_submodel_info.microenvironment_variables.push_back( "danger signals" ); 
	//epithelium_submodel_info.microenvironment_variables.push_back( "interferon 1" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "pro-inflammatory cytokine" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "chemokine" ); 
	epithelium_submodel_info.microenvironment_variables.push_back( "anti-inflammatory cytokine" );
	//epithelium_submodel_info.microenvironment_variables.push_back( "pro-pyroptosis cytokine" );
		// what custom data do I need? 
	//epithelium_submodel_info.cell_variables.push_back( "something" ); 
		// register the submodel  
	epithelium_submodel_info.register_model();	
		// set functions for the corresponding cell definition 
	pCD = find_cell_definition( "skin cell" ); 
	pCD->functions.update_phenotype = epithelium_submodel_info.phenotype_function;
	pCD->functions.custom_cell_rule = epithelium_submodel_info.mechanics_function;
	pCD->functions.contact_function = epithelium_contact_function; 
	
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


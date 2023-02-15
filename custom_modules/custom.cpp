/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                     #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
// vector of valid positions
std::vector<std::vector <double>> valid_position;
// Cells that will receive a nudge to come back to the domain
std::vector<Cell*> cells_to_move_from_edge;
extern LymphNode lympNode;
extern std::vector<int> vascularized_voxel_indices;

void create_cell_types( void )
{
	// set the random seed
	SeedRandom( parameters.ints("random_seed") );

	/*
	Put any modifications to default cell definition here if you
	want to have "inherited" by other cell types.

	This is a good place to set default functions.
	*/
	initialize_default_cell_definition();
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL;
	cell_defaults.functions.custom_cell_rule = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*
	This parses the cell definitions in the XML config file.
	*/
	// Variables of plastoelastic mechanics
	std::vector<double> myvec(3,0.0);
	cell_defaults.custom_data.add_vector_variable( "ECM_attachment_point", "micron", myvec );
	cell_defaults.custom_data.add_vector_variable( "mechanical_strain_displacement" , "micron", myvec );
	// assume elastic movement on the order of 10 min at maximum 10 micron elongation
	cell_defaults.custom_data.add_variable( "spring_constant", "1/min" , 0.05 );  // 0.05    // (1.0/10.0) * (1.0/10.0)
	// assume plastic movement on the order of 1 day at maximum 10 micron elongation
	cell_defaults.custom_data.add_variable( "mechanical_relaxation_rate", "1/min" , 0.0005 );  // 0.0005  // (1.0/10.0) * (1.0/(24.0*60.0)
	cell_defaults.custom_data.add_variable( "mechanical_strain", "micron" , 0.0 );

	initialize_cell_definitions_from_pugixml();

	/*
	Put any modifications to individual cell definitions here.

	This is a good place to set custom functions.
	*/

	// register the submodels
	// (which ensures that the cells have all the internal variables they need)
	immune_submodels_setup();
	epithelium_submodel_setup();
	cancer_submodel_setup();

	submodel_registry.display( std::cout );

	/*
	This builds the map of cell definitions and summarizes the setup.
	*/
	build_cell_definitions_maps();
	display_cell_definitions( std::cout );

	return;
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl;
		default_microenvironment_options.simulate_2D = true;
	}

	initialize_microenvironment();

	return;
}

void setup_tissue( void )
{
	choose_initialized_voxels();
	// create some cells near the origin

	Cell* pC;

	// hexagonal cell packing
	Cell_Definition* pCD = find_cell_definition("lung cell");

	double cell_radius = pCD->phenotype.geometry.radius;
	double spacing = 0.95 * cell_radius * 2.0;

	double x_min = microenvironment.mesh.bounding_box[0] + cell_radius;
	double x_max = microenvironment.mesh.bounding_box[3] - cell_radius;

	double y_min = microenvironment.mesh.bounding_box[1] + cell_radius;
	double y_max = microenvironment.mesh.bounding_box[4] - cell_radius;

	double x = x_min;
	double y = y_min;

	double center_x = 0.5*( x_min + x_max );
	double center_y = 0.5*( y_min + y_max );

	double triangle_stagger = sqrt(3.0) * spacing * 0.5;

	// find the cell nearest to the center
	double nearest_distance_squared = 9e99;
	Cell* pNearestCell = NULL;

	std::vector<double> temp_position(3, 0.0);

	int n = 0;
	while( y < y_max )
	{
		while( x < x_max )
		{
			temp_position[0] = x; temp_position[1] = y; temp_position[2] = 0.0;
			valid_position.push_back(temp_position);

			double dx = x - center_x;
			double dy = y - center_y;

			double temp = dx*dx + dy*dy;
			if( temp < nearest_distance_squared )
			{
				nearest_distance_squared = temp;
				pNearestCell = pC;
			}
			x += spacing;
		}
		x = x_min;

		n++;
		y += triangle_stagger;
		// in odd rows, shift
		if( n % 2 == 1 )
		{
			x += 0.5 * spacing;
		}
	}
	int Max_number_of_cell = valid_position.size();
	// extern double GridCOUNT;
  lympNode.GridCOUNT = Max_number_of_cell;
	// place immune cells
	initial_immune_cell_placement();

	static int ECM_attachment_point_index = pCD->custom_data.find_vector_variable_index( "ECM_attachment_point" );
	// cancer cells
	if ( parameters.doubles("time_add_cancer_cell") == 0.0 )
	{
		for ( int i=0; i < parameters.ints("number_of_cancer_cells") ;i++){
			pC = create_cell( get_cell_definition("cancer cell" ) );
			int index_sample = (int) ( UniformRandom() * (valid_position.size()-1) );
			pC->assign_position( valid_position[index_sample] );
			pC->custom_data.vector_variables[ECM_attachment_point_index].value = pC->position;
			valid_position.erase (valid_position.begin()+index_sample);
			//std::cout << "SIZE: " << valid_position.size() << " Index: " << index_sample <<  std::endl;
		}
	}
	// Ephitelium cells
  	for ( int i=0; i < parameters.doubles("cell_confluence_lung_cells")*Max_number_of_cell;i++){
		pC = create_cell( get_cell_definition("lung cell" ) );
		int index_sample = (int) ( UniformRandom() * (valid_position.size()-1) );
		pC->assign_position( valid_position[index_sample] );
		pC->custom_data.vector_variables[ECM_attachment_point_index].value = pC->position;
		valid_position.erase (valid_position.begin()+index_sample);
		if (valid_position.size() == 0) break;
	}

	return;
}

std::vector<std::string> epithelium_coloring_function( Cell* pCell )
{
	std::vector<std::string> output( 4, "black" );

	static int color_index =
	cell_defaults.custom_data.find_variable_index( parameters.strings["color_variable"].value );

	// color by assembled generation
	if( pCell->phenotype.death.dead == false )
	{
		// find fraction of max generation load
		double v = pCell->custom_data[ color_index ] ;

		double interpolation = 0;
		if( v < 10.0 )
		{ interpolation = 1.0 - v*0.05; }
		else
		{ interpolation = 0.5; }

		int red = 0;
		int green = 0;
		int blue = (int) floor( 255.0 * interpolation );

		char color [1024];
		sprintf( color, "rgb(%u,%u,%u)" , red,green,blue );

		output[0] = color;
		output[2] = color;
		output[3] = color;
	}

	return output;
}

std::vector<std::string> tissue_coloring_function( Cell* pCell )
{
	static int lung_epithelial_type = get_cell_definition( "lung cell" ).type;
	static int cancer_type = get_cell_definition( "cancer cell" ).type;

	static int CD8_Tcell_type = get_cell_definition( "CD8 Tcell" ).type;
	static int Macrophage_type = get_cell_definition( "macrophage" ).type;
	static int DC_type = get_cell_definition( "DC" ).type;
	static int CD4_Tcell_type = get_cell_definition( "CD4 Tcell" ).type;

	// start with white

	std::vector<std::string> output = {"white", "black", "white" , "white" };
	// false_cell_coloring_cytometry(pCell);

	if( pCell->phenotype.death.dead == true )
	{
		if( pCell->type != lung_epithelial_type && pCell->type != cancer_type )
		{
			output[0] = parameters.strings("apoptotic_immune_color");
			output[2] = output[0];
			output[3] = output[0];
			return output;
		}
		else{
			if( pCell->type == lung_epithelial_type)
			{
				output[0] = parameters.strings("apoptotic_epithelium_color");
				output[2] = output[0];
				output[3] = output[0];
				return output;
			}
			else{
				output[0] = parameters.strings("apoptotic_cancer_color");
				output[2] = output[0];
				output[3] = output[0];
			}
		}
	}

	if( pCell->phenotype.death.dead == false && pCell->type == lung_epithelial_type )
	{
		output[0] = "rgb(0,0,255)";
		output[2] = "rgb(0,0,255)";
		output[3] = "rgb(0,0,255)";
		// output = epithelium_coloring_function(pCell);
		return output;
	}

	if( pCell->phenotype.death.dead == false && pCell->type == cancer_type )
	{
		output[0] = parameters.strings("cancer_color");
		output[2] = output[0];
		output[3] = output[0];
		return output;
	}

	if( pCell->phenotype.death.dead == false && pCell->type == CD8_Tcell_type )
	{
		output[0] = parameters.strings("CD8_Tcell_color");
		output[2] = output[0];
		output[3] = output[0];
		return output;
	}

	// (Adrianne) adding CD4 T cell colouring
	if( pCell->phenotype.death.dead == false && pCell->type == CD4_Tcell_type )
	{
		output[0] = parameters.strings("CD4_Tcell_color");
		output[2] = output[0];
		output[3] = output[0];
		return output;
	}

	if( pCell->phenotype.death.dead == false && pCell->type == Macrophage_type )
	{
		std::string color = parameters.strings("Macrophage_color");
		if( pCell->custom_data["activated_immune_cell" ] > 0.5 )
		{ color = parameters.strings("activated_macrophage_color"); }

		if( pCell->phenotype.volume.total> pCell->custom_data["threshold_macrophage_volume"] )// macrophage exhausted
		{ color = parameters.strings("exhausted_macrophage_color"); }
		else if( pCell->custom_data["ability_to_phagocytose_cancer_cell"] == 1)// macrophage has been activated to kill cancer cells by T cell
		{ color = parameters.strings("hyperactivated_macrophage_color"); }

		output[0] = color;
		output[2] = output[0];
		output[3] = output[0];
		return output;
	}

	//(Adrianne) adding colour for DCs
	if( pCell->phenotype.death.dead == false && pCell->type == DC_type )
	{
		std::string color = parameters.strings("DC_color");
		if( pCell->custom_data["activated_immune_cell" ] > 0.5 )
		{ color = parameters.strings("activated_DC_color"); }

		output[0] = color;
		output[2] = output[0];
		output[3] = output[0];
		return output;
	}

	return output;
}

bool Write_SVG_circle_opacity( std::ostream& os, double center_x, double center_y, double radius, double stroke_size,	std::string stroke_color , std::string fill_color , double opacity )
{
	os << "  <circle cx=\"" << center_x << "\" cy=\"" << center_y << "\" r=\"" << radius << "\" stroke-width=\"" << stroke_size
	<< "\" stroke=\"" << stroke_color << "\" fill=\"" << fill_color
	<< "\" fill-opacity=\"" << opacity << "\"/>" << std::endl;
	return true;
}

void SVG_plot_custom( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*) )
{
		static double X_lower = M.mesh.bounding_box[0];
		static double X_upper = M.mesh.bounding_box[3];

		static double Y_lower = M.mesh.bounding_box[1];
		static double Y_upper = M.mesh.bounding_box[4];

		static double plot_width = X_upper - X_lower;
		static double plot_height = Y_upper - Y_lower;

		static double font_size = 0.025 * plot_height; // PhysiCell_SVG_options.font_size;
		static double top_margin = font_size*(.2+1+.2+.9+.5 );

		static double epithelial_opacity = parameters.doubles("epithelial_opacity");
		static double non_epithelial_opacity = parameters.doubles("non_epithelial_opacity");

		// open the file, write a basic "header"
		std::ofstream os( filename , std::ios::out );
		if( os.fail() )
		{
				std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl;

				std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
				<< "Check to make sure your save directory exists. " << std::endl << std::endl
				<< "I'm going to exit with a crash code of -1 now until " << std::endl
				<< "you fix your directory. Sorry!" << std::endl << std::endl;
				exit(-1);
		}

		Write_SVG_start( os, plot_width , plot_height + top_margin );

		// draw the background
		Write_SVG_rect( os , 0 , 0 , plot_width, plot_height + top_margin , 0.002 * plot_height , "white", "white" );

		// write the simulation time to the top of the plot

		char* szString;
		szString = new char [1024];

		int total_cell_count = all_cells->size();

		double temp_time = time;

		std::string time_label = formatted_minutes_to_DDHHMM( temp_time );

		sprintf( szString , "Current time: %s, z = %3.2f %s", time_label.c_str(),
		z_slice , PhysiCell_SVG_options.simulation_space_units.c_str() );
		Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1),
		font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
		sprintf( szString , "%u agents" , total_cell_count );
		Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1+.2+.9),
		0.95*font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );

		delete [] szString;


		// add an outer "g" for coordinate transforms

		os << " <g id=\"tissue\" " << std::endl
		<< "    transform=\"translate(0," << plot_height+top_margin << ") scale(1,-1)\">" << std::endl;

		// prepare to do mesh-based plot (later)

		double dx_stroma = M.mesh.dx;
		double dy_stroma = M.mesh.dy;

		os << "  <g id=\"ECM\">" << std::endl;

		int ratio = 1;
		double voxel_size = dx_stroma / (double) ratio ;

		double half_voxel_size = voxel_size / 2.0;
		double normalizer = 78.539816339744831 / (voxel_size*voxel_size*voxel_size);

		os << "  </g>" << std::endl;

		static Cell_Definition* pEpithelial = find_cell_definition( "lung cell" );
		static Cell_Definition* pcancer = find_cell_definition( "cancer cell" );

		// plot intersecting epithelial and cancer cells
		os << "  <g id=\"cells\">" << std::endl;
		for( int i=0 ; i < total_cell_count ; i++ )
		{
				Cell* pC = (*all_cells)[i]; // global_cell_list[i];

				static std::vector<std::string> Colors;
				if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius
				&& (pC->type == pEpithelial->type || pC->type == pcancer->type) )
				{
						double r = pC->phenotype.geometry.radius ;
						double rn = pC->phenotype.geometry.nuclear_radius ;
						double z = fabs( (pC->position)[2] - z_slice) ;

						Colors = cell_coloring_function( pC );

						os << "   <g id=\"cell" << pC->ID << "\">" << std::endl;

						// figure out how much of the cell intersects with z = 0

						double plot_radius = sqrt( r*r - z*z );

						Write_SVG_circle_opacity( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower,
						plot_radius , 0.5, Colors[1], Colors[0] , epithelial_opacity );
						os << "   </g>" << std::endl;
				}
		}

		//plot intersecting non=epithelial and non-cancer cells
		for( int i=0 ; i < total_cell_count ; i++ )
		{
				Cell* pC = (*all_cells)[i]; // global_cell_list[i];

				static std::vector<std::string> Colors;
				if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius
				&& pC->type != pEpithelial->type && pC->type != pcancer->type)
				{
						double r = pC->phenotype.geometry.radius ;
						double rn = pC->phenotype.geometry.nuclear_radius ;
						double z = fabs( (pC->position)[2] - z_slice) ;

						Colors = cell_coloring_function( pC );

						os << "   <g id=\"cell" << pC->ID << "\">" << std::endl;

						// figure out how much of the cell intersects with z = 0

						double plot_radius = sqrt( r*r - z*z );

						Write_SVG_circle_opacity( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower,
						plot_radius , 0.5, Colors[1], Colors[0] , non_epithelial_opacity );
						os << "   </g>" << std::endl;
				}
		}

		os << "  </g>" << std::endl;

		// end of the <g ID="tissue">
		os << " </g>" << std::endl;

		// draw a scale bar

		double bar_margin = 0.025 * plot_height;
		double bar_height = 0.01 * plot_height;
		double bar_width = PhysiCell_SVG_options.length_bar;
		double bar_stroke_width = 0.001 * plot_height;

		std::string bar_units = PhysiCell_SVG_options.simulation_space_units;
		// convert from micron to mm
		double temp = bar_width;

		if( temp > 999 && std::strstr( bar_units.c_str() , PhysiCell_SVG_options.mu.c_str() )   )
		{
				temp /= 1000;
				bar_units = "mm";
		}
		// convert from mm to cm
		if( temp > 9 && std::strcmp( bar_units.c_str() , "mm" ) == 0 )
		{
				temp /= 10;
				bar_units = "cm";
		}

		szString = new char [1024];
		sprintf( szString , "%u %s" , (int) round( temp ) , bar_units.c_str() );

		Write_SVG_rect( os , plot_width - bar_margin - bar_width  , plot_height + top_margin - bar_margin - bar_height ,
		bar_width , bar_height , 0.002 * plot_height , "rgb(255,255,255)", "rgb(0,0,0)" );
		Write_SVG_text( os, szString , plot_width - bar_margin - bar_width + 0.25*font_size ,
		plot_height + top_margin - bar_margin - bar_height - 0.25*font_size ,
		font_size , PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );

		delete [] szString;

		// plot runtime
		szString = new char [1024];
		RUNTIME_TOC();
		std::string formatted_stopwatch_value = format_stopwatch_value( runtime_stopwatch_value() );
		Write_SVG_text( os, formatted_stopwatch_value.c_str() , bar_margin , top_margin + plot_height - bar_margin , 0.75 * font_size ,
		PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
		delete [] szString;

		// draw a box around the plot window
		Write_SVG_rect( os , 0 , top_margin, plot_width, plot_height , 0.002 * plot_height , "rgb(0,0,0)", "none" );

		// close the svg tag, close the file
		Write_SVG_end( os );
		os.close();

		return;
}

void divide_custom_data()
{
	static int cancer_type = get_cell_definition( "cancer cell" ).type;
	static double tolerance = 0.001;

  #pragma omp parallel for
	for( int i=0; i < (*all_cells).size() ;i++ )
	{
		Cell* pCell;
		pCell = (*all_cells)[i];

		// if cell is dead or is not cancer cell, skip it
		if( pCell->phenotype.death.dead == true || pCell->type != cancer_type)
	    { continue; }

		static int last_cycle_index = pCell->custom_data.find_variable_index( "last_cycle_entry_time");
		static int generation_index = pCell->custom_data.find_variable_index( "division_generation" );
		if( pCell->phenotype.cycle.data.elapsed_time_in_phase < tolerance &&
			fabs( PhysiCell_globals.current_time - pCell->custom_data[ last_cycle_index ] ) >= phenotype_dt )
		{
			// add generation by 1
			pCell->custom_data[ generation_index ] += 1;
			pCell->custom_data[ last_cycle_index ] = PhysiCell_globals.current_time;
		}
	}

	return;
}

std::vector<double> set_nudge_from_edge( Cell* pC , double tolerance ) // return {push_x,push_y,push_z} of direction to nudge cell
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
	static int cancer_type = get_cell_definition( "cancer cell" ).type;
	for( int n=0 ; n < cells_to_move_from_edge.size(); n++ )
	{
		Cell* pC = cells_to_move_from_edge[n];
		if ( pC->type == cancer_type){
			delete_cell(pC);
		}else{
			nudge_out_of_bounds_cell( pC , 10.0 );
		}
	}
	return;
}

void clear_cells_to_move_from_edge( void )
{
	cells_to_move_from_edge.clear();
}

void include_tumor_cells(void)
{
	if (  PhysiCell_globals.current_time == 0 || fabs( PhysiCell_globals.current_time - parameters.doubles("time_add_cancer_cell")  ) >= 0.01 * diffusion_dt )
		return;
	Cell* pC;
	Cell_Definition* pCD = find_cell_definition( "cancer cell" );
	static int ECM_attachment_point_index = pCD->custom_data.find_vector_variable_index( "ECM_attachment_point" );
	for ( int i=0; i < parameters.ints("number_of_cancer_cells") ;i++){
		pC = create_cell( get_cell_definition("cancer cell" ) );
		int index_sample = (int) ( UniformRandom() * (valid_position.size()-1) );
		pC->assign_position( valid_position[index_sample] );
		pC->custom_data.vector_variables[ECM_attachment_point_index].value = pC->position;
		valid_position.erase (valid_position.begin()+index_sample);
	}
}

void vaccine(void)
{
	if ( parameters.ints("number_of_cells_vaccine") == 0 ) return; // No vaccine
	if (  fabs(PhysiCell_globals.current_time - (parameters.doubles("time_vaccine") + parameters.ints("count_vaccine")*parameters.doubles("vaccine_interval"))) >= 0.01 * diffusion_dt ) // Check right time
		return;
	parameters.ints("count_vaccine")++;
	std::cout<< "Take shot: " << PhysiCell_globals.current_time << " min - dose: " << parameters.ints("number_of_cells_vaccine") << " cells\n";

	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Xmax = microenvironment.mesh.bounding_box[3];
	double Ymax = microenvironment.mesh.bounding_box[4];
	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	Cell* pC;
	Cell_Definition* pCD = find_cell_definition( "cancer cell" );

	for ( int i=0; i < parameters.ints("number_of_cells_vaccine") ;i++){
		std::vector<double> position = {0,0,0};
		position[0] = Xmin + UniformRandom()*Xrange;
		position[1] = Ymin + UniformRandom()*Yrange;
		pC = create_cell( get_cell_definition("cancer cell" ) );
		pC->assign_position( position );
		//pC->custom_data.vector_variables[ECM_attachment_point_index].value = pC->position;
		static int apoptosis_index = 	pC->phenotype.death.find_death_model_index( "Apoptosis" );
		static int debris_index = microenvironment.find_density_index( "debris" );
		pC->start_death( apoptosis_index );
		pC->phenotype.volume.total = 478;
		pC->phenotype.volume.nuclear = 47.8;
		pC->phenotype.geometry.update( pC, pC->phenotype, 0.0 );
		pC->phenotype.cycle.data.transition_rate( 0, 1 ) = 0.0;
		pC->phenotype.secretion.secretion_rates[debris_index] = pC->custom_data["debris_secretion_rate"];
		pC->functions.volume_update_function = NULL;
		pC->functions.update_phenotype = NULL;
	}
}

std::vector<int> cell_count( void )
{
	static int lung_cell_type = get_cell_definition( "lung cell" ).type;
	static int cancer_cell_type = get_cell_definition( "cancer cell" ).type;
	static int CD8_type = get_cell_definition( "CD8 Tcell" ).type;
	static int CD4_type = get_cell_definition( "CD4 Tcell" ).type;
	static int macrophage_type = get_cell_definition( "macrophage" ).type;
	static int DC_type = get_cell_definition( "DC" ).type;

	std::vector<int> NumberofCells(16,0);
	for (int i=0; i < (*all_cells).size(); i++)
	{
		if( (*all_cells)[i]->phenotype.death.dead == true )
		{
			if( (*all_cells)[i]->type == lung_cell_type ) NumberofCells[1]++; // Dead lung
			else if( (*all_cells)[i]->type == cancer_cell_type ) NumberofCells[2]++; // Dead cancer
			else if( (*all_cells)[i]->type == DC_type ) NumberofCells[3]++; // Dead DC
			else if( (*all_cells)[i]->type == macrophage_type ) NumberofCells[4]++; // Dead macrophage
			else if( (*all_cells)[i]->type == CD4_type ) NumberofCells[5]++; // Dead CD4
			else if( (*all_cells)[i]->type == CD8_type ) NumberofCells[6]++; // Dead CD8
		}
		else if( (*all_cells)[i]->type == lung_cell_type )
		{
			NumberofCells[0]++;
		}
		else if ( (*all_cells)[i]->type == cancer_cell_type )
		{
			NumberofCells[7]++;
		}
		else if ( (*all_cells)[i]->type == CD8_type )
		{
			NumberofCells[15]++;
		}
		else if ( (*all_cells)[i]->type == macrophage_type )
		{
			if ( (*all_cells)[i]->custom_data["ability_to_phagocytose_cancer_cell"] == 1 )
			{	NumberofCells[11]++; } //hyperactivated macrophage
			else
			{
				if ( (*all_cells)[i]->phenotype.volume.total> (*all_cells)[i]->custom_data["threshold_macrophage_volume"] ){ NumberofCells[10]++; } //exhausted macrophage
				else
				{
					if ( (*all_cells)[i]->custom_data["activated_immune_cell" ] > 0.5 ){ NumberofCells[14]++; } //activated macrophage
					else NumberofCells[9]++; //inactivated macrophage
				}
			}
		}
		else if ( (*all_cells)[i]->type == DC_type )
		{
			if ( (*all_cells)[i]->custom_data["activated_immune_cell" ] > 0.5 ){ NumberofCells[13]++; } //activated DC
			else NumberofCells[8]++; //inactivated DC
		}
		else if ( (*all_cells)[i]->type == CD4_type )
		{
			NumberofCells[12]++;
		}
	}
	return NumberofCells;
}

void Binary2Color(const std::vector<bool> &BinaryVector, int &Rvalue, int &Gvalue, int &Bvalue)
{
  // -------- Red channel -------------
  // 0 - live lung cell;              1 - dead lung cell;
  // 2 - dead cancer cell;          3 - dead DC
  // 4 - dead macrophage;             5 - dead CD4 T cell;
  // 6 - dead CD8 T cell;             7 - live cancer cell;
  // -------- Green channel -------------
  // 8 - inactivated DC;              9 - inactivated macrophage;
  // 10 - exhausted macrophage;       11 - hyperactivated macrophage;
  // 12 - live CD4 T cell;            13 - activated DC;
  // 14 - activated macrophage;       15 - live CD8 T cell;
  // -------- Blue channel -------------
  // 16 - debris >= 0.00 and < 0.25;  17 - debris >= 0.25 and < 0.50;
  // 18 - debris >= 0.50 and < 0.75;  19 - debris >= 0.75;
  // 20 - TNF >= 0.00 and < 0.25;     21 - TNF >= 0.25 and < 0.50;
  // 22 - TNF >= 0.50 and < 0.75;     23 - TNF >= 0.75;
  Rvalue = 0; Gvalue = 0; Bvalue = 0;
  // Red channel
  for (unsigned int i = 0; i < 8; i++){
      Rvalue += BinaryVector[i]*pow(2,i);
  }
  // Green channel
  for (unsigned int i = 8; i < 16; i++){
      Gvalue += BinaryVector[i]*pow(2,(i-8));
  }
  // Blue channel
  for (unsigned int i = 16; i < 24; i++){
      Bvalue += BinaryVector[i]*pow(2,(i-16));
  }
}

void SetBinaryVector(const std::vector<double> &Pos, const int type, std::vector<std::vector<binaryVec>> &PixelsBinary, BMP &Image )
{
	const double x_min = microenvironment.mesh.bounding_box[0];
	const double x_max = microenvironment.mesh.bounding_box[3];
	const double y_min = microenvironment.mesh.bounding_box[1];
	const double y_max = microenvironment.mesh.bounding_box[4];

  const double deltaX = (x_max - x_min)/Image.TellWidth();
  const double deltaY = (y_max - y_min)/Image.TellHeight();
  // Check cells inside of pixels - type is the binary position
  for(unsigned int i = 0; i < Image.TellWidth(); i++){
    for (unsigned int j = 0; j < Image.TellHeight(); j++){
        if ( Pos[0] >= x_min + deltaX*i && Pos[0] < x_min + deltaX*(i+1) &&
        Pos[1] <= y_max - deltaY*j && Pos[1] > y_max - deltaY*(j+1) )
        {
          PixelsBinary[i][j].binaryVector[type] = true;
          //cout <<"(" << Pos[0] << ", " << Pos[1] << ") - (" << i << ", " << j <<") - "<< type << endl;
          return;
        }

    }
  }
}

void SetBinaryVectorSubs(const std::vector<double> &Pos, const double debris,const double TNF, std::vector<std::vector<binaryVec>> &PixelsBinary, BMP &Image )
{
  const double x_min = microenvironment.mesh.bounding_box[0];
	const double x_max = microenvironment.mesh.bounding_box[3];
	const double y_min = microenvironment.mesh.bounding_box[1];
	const double y_max = microenvironment.mesh.bounding_box[4];

  const double deltaX = (x_max - x_min)/Image.TellWidth();
  const double deltaY = (y_max - y_min)/Image.TellHeight();
  // Check center of voxels inside of pixels
  for(unsigned int i = 0; i < Image.TellWidth(); i++){
    for (unsigned int j = 0; j < Image.TellHeight(); j++){
        if ( Pos[0] >= x_min + deltaX*i && Pos[0] < x_min + deltaX*(i+1) &&
        Pos[1] <= y_max - deltaY*j && Pos[1] > y_max - deltaY*(j+1) )
        {
          if (debris >= 0.00 &&  debris < 0.25) PixelsBinary[i][j].binaryVector[16] = true;
          if (debris >= 0.25 &&  debris < 0.50) PixelsBinary[i][j].binaryVector[17] = true;
          if (debris >= 0.50 &&  debris < 0.75) PixelsBinary[i][j].binaryVector[18] = true;
          if (debris >= 0.75 ) PixelsBinary[i][j].binaryVector[19] = true;
          if (TNF >= 0.00 &&  TNF < 0.25) PixelsBinary[i][j].binaryVector[20] = true;
          if (TNF >= 0.25 &&  TNF < 0.50) PixelsBinary[i][j].binaryVector[21] = true;
          if (TNF >= 0.50 &&  TNF < 0.75) PixelsBinary[i][j].binaryVector[22] = true;
          if (TNF >= 0.75 ) PixelsBinary[i][j].binaryVector[23] = true;

          //cout <<"(" << Pos[0] << ", " << Pos[1] << ") - (" << i << ", " << j <<") - "<< type << endl;
          return;
        }

    }
  }
}

void GenerateBitmap(const char* filename)
{
	BMP Image;
  // Set size to width x height
  Image.SetSize(parameters.ints("bitmap_width"),parameters.ints("bitmap_height"));
  // Set its color depth to 24-bits
  Image.SetBitDepth(24);
  std::vector<std::vector<binaryVec>> PixelsBinary(Image.TellHeight(), std::vector<binaryVec>(Image.TellWidth()));

	static int lung_type = get_cell_definition( "lung cell" ).type;
	static int cancer_type = get_cell_definition( "cancer cell" ).type;
	static int CD8_type = get_cell_definition( "CD8 Tcell" ).type;
	static int CD4_type = get_cell_definition( "CD4 Tcell" ).type;
	static int macrophage_type = get_cell_definition( "macrophage" ).type;
	static int DC_type = get_cell_definition( "DC" ).type;
	for (int i=0; i < (*all_cells).size(); i++)
	{
		if( (*all_cells)[i]->type == lung_type ){
			if ( (*all_cells)[i]->phenotype.death.dead == false ) SetBinaryVector((*all_cells)[i]->position,0,PixelsBinary,Image); //bit 0 - live lung cell
			else SetBinaryVector((*all_cells)[i]->position,1,PixelsBinary,Image); //bit 1 - dead lung cell
		}
		else if ( (*all_cells)[i]->type == cancer_type ){
			if ( (*all_cells)[i]->phenotype.death.dead == true ) SetBinaryVector((*all_cells)[i]->position,2,PixelsBinary,Image); //bit 2 - dead cancer cell
			else SetBinaryVector((*all_cells)[i]->position,7,PixelsBinary,Image); //bit 7 - live cancer cell
		}
		else if ( (*all_cells)[i]->type == DC_type ){
			if ( (*all_cells)[i]->phenotype.death.dead == true ) SetBinaryVector((*all_cells)[i]->position,3,PixelsBinary,Image); //bit 3 - dead DC
			else if ( (*all_cells)[i]->custom_data["activated_immune_cell" ] > 0.5 ){
				SetBinaryVector((*all_cells)[i]->position,13,PixelsBinary,Image); //bit 13 - activated DC
			}else SetBinaryVector((*all_cells)[i]->position,8,PixelsBinary,Image); //bit 8 - inactivated DC
		}
		else if ( (*all_cells)[i]->type == macrophage_type ){
			if ( (*all_cells)[i]->phenotype.death.dead == true ) SetBinaryVector((*all_cells)[i]->position,4,PixelsBinary,Image); //bit 4 - dead macrophage
			else if ( (*all_cells)[i]->custom_data["ability_to_phagocytose_cancer_cell"] == 1 ){
				SetBinaryVector((*all_cells)[i]->position,11,PixelsBinary,Image); //bit 11 - hyperactivated macrophage
			}else if ( (*all_cells)[i]->phenotype.volume.total> (*all_cells)[i]->custom_data["threshold_macrophage_volume"] ){
				SetBinaryVector((*all_cells)[i]->position,10,PixelsBinary,Image); //bit 10 - exhausted macrophage
			}else if ( (*all_cells)[i]->custom_data["activated_immune_cell" ] > 0.5 ){
				SetBinaryVector((*all_cells)[i]->position,14,PixelsBinary,Image); //bit 14 - activated macrophage
			}else SetBinaryVector((*all_cells)[i]->position,9,PixelsBinary,Image); //bit 9 - inactivated macrophage
		}
		else if( (*all_cells)[i]->type == CD4_type ){
			if ( (*all_cells)[i]->phenotype.death.dead == true ) SetBinaryVector((*all_cells)[i]->position,5,PixelsBinary,Image); //bit 5 - dead CD4 T cell
			else SetBinaryVector((*all_cells)[i]->position,12,PixelsBinary,Image); //bit 12 - live CD4 T cell
		}
		else if( (*all_cells)[i]->type == CD8_type ){
			if ( (*all_cells)[i]->phenotype.death.dead == true ) SetBinaryVector((*all_cells)[i]->position,6,PixelsBinary,Image); //bit 6 - dead CD8 T cell
			else SetBinaryVector((*all_cells)[i]->position,15,PixelsBinary,Image); //bit 15 - live CD8 T cell
		}
	}

	// Substrate in blue channel
	static int TNF_index = microenvironment.find_density_index( "TNF");
	static int debris_index = microenvironment.find_density_index( "debris");
	for( unsigned int n=0; n < microenvironment.number_of_voxels() ; n++ ){
		SetBinaryVectorSubs(microenvironment.mesh.voxels[n].center,microenvironment.nearest_density_vector(n)[debris_index],microenvironment.nearest_density_vector(n)[TNF_index],PixelsBinary,Image);
		//std::cout << microenvironment.mesh.voxels[n].center <<" "<< microenvironment.nearest_density_vector(n)[debris_index] <<" "<< microenvironment.nearest_density_vector(n)[TNF_index] << std::endl;
	}

	// Write bitmap file
	int Rvalue, Gvalue, Bvalue;
  for(unsigned int i = 0; i < Image.TellWidth(); i++)
    for (unsigned int j = 0; j < Image.TellHeight(); j++){
      Binary2Color(PixelsBinary[i][j].binaryVector,Rvalue,Gvalue,Bvalue);
      Image(i,j)->Red = Rvalue;
      Image(i,j)->Green = Gvalue;
      Image(i,j)->Blue = Bvalue;
      // std::cout << "Pixels - (" << i << ", " << j << ") - binary: [ ";
      // for (unsigned int k = 0; k < PixelsBinary[i][j].binaryVector.size(); k++)
      //   std::cout << PixelsBinary[i][j].binaryVector[k] << " ";
      // std::cout << " ] - R: " << Rvalue << " G: " << Gvalue << " B: " << Bvalue << std::endl;
    }

	// Write bitmap file
	Image.WriteToFile( filename );


}

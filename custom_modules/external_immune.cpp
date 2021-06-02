#include "./external_immune.h"

using namespace PhysiCell;

std::string external_immune_version = "0.4.0";

Submodel_Information external_immune_info;

// generator sample from distribution
extern std::default_random_engine generator;

extern AntigenLibrary AntigenLib;

void external_immune_model_setup( void )
{
	// set version
	external_immune_info.name = "external immune";
	external_immune_info.version = external_immune_version;
	// set functions
	external_immune_info.main_function = external_immune_model;
	external_immune_info.phenotype_function = NULL;
	external_immune_info.mechanics_function = NULL;
	// what microenvironment variables do I need?

	// what custom data do I need?
	//external_immune_info.parameters.doubles.push_back( "DM" );
	//external_immune_info.parameters.doubles.push_back( "TC" );

	// submodel_registry.register_model( internal_pathogen_dynamics_info );
	external_immune_info.register_model();

	return;
}

void external_immune_model( double dt )
{
	AntigenLib.update_collection(PhysiCell_globals.current_time, 720.0); // update collection of antigens
	// bookkeeping -- find microenvironment variables we need

	extern double DM;
	extern double TC;
	extern double TCt;
	extern double GridCOUNT;
	static double dC = parameters.doubles( "TC_death_rate" );
	static double pT1 = parameters.doubles( "max_activation_TC" );
	static double pT2 = parameters.doubles( "half_max_activation_TC" );
	static double dT1 = parameters.doubles( "max_clearance_TC" );
	static double dT2 = parameters.doubles( "half_max_clearance_TC" );
	static double Tc0 = parameters.doubles( "TC_population_threshold" );
	static double immunevolume = 1;
	static double dDm = parameters.doubles( "DM_decay" );
	static double ro = 1;
	static double CD8_Tcell_recruitment_rate = parameters.doubles( "T_Cell_Recruitment" );

	//double lypmh_scale = GridCOUNT / 3e6; // default 5e6
	double lypmh_scale = 0.02; // old version

	// actual model goes here

	double x[4][3]={{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};//initialize x
	double f[4][3]={{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};//initialize f
	int j;

	// TC update
	double dR_TC = dC * Tc0;
	x[0][0] = DM/lypmh_scale;
	x[0][1] = TC; //initial values
	x[0][2] = TCt/lypmh_scale;

	for(j = 0; j < 4; j++){
		f[j][0] = {-dDm*x[j][0]}; //define function
		f[j][1] = {dR_TC-dC*x[j][1]+pT1*x[j][0]*x[j][1]/(x[j][0]+pT2)-dT1*x[j][0]*x[j][1]/(x[j][0]+dT2)};// /* TEST */ + (1e6 - x[j][1])/1e6};
		f[j][2] = {CD8_Tcell_recruitment_rate*x[j][1]}; //define function

		if (j== 0 || j==1){
			x[j+1][0]=x[0][0]+dt/2*f[j][0]; //first and second x approximations
			x[j+1][1]=x[0][1]+dt/2*f[j][1]; //first and second x approximations
			x[j+1][2]=x[0][2]+dt/2*f[j][2]; //first and second x approximations
		}
		if (j== 2){
			x[j+1][0]=x[0][0]+dt*f[j][0]; //third approximation
			x[j+1][1]=x[0][1]+dt*f[j][1]; //third approximation
			x[j+1][2]=x[0][2]+dt*f[j][2]; //third approximation
		}
	}

	//std::cout << dt*(f[0][0]/6+f[1][0]/3+f[2][0]/3+f[3][0]/6) << std::endl;

	DM=(x[0][0]+dt*(f[0][0]/6+f[1][0]/3+f[2][0]/3+f[3][0]/6))*lypmh_scale;
	TC=x[0][1]+dt*(f[0][1]/6+f[1][1]/3+f[2][1]/3+f[3][1]/6);
	TCt=(x[0][2]+dt*(f[0][2]/6+f[1][2]/3+f[2][2]/3+f[3][2]/6))*lypmh_scale;


	return;
}

// Methods from neoantigens_library

bool AntigenLibrary::check_library(const std::vector<double> &antigen_signature)
{
	// Check collection
	for( int i=0; i < Collection.size() ;i++ )
	{
		for( int j=0; j < Collection[i].signature.size() ;j++ )
			if ( Collection[i].signature[j] != antigen_signature[j] ) break;
			else if( j == Collection[i].signature.size()-1 ){ Collection[i].probability_weight++; return true;}
	}
	// Check temporary collection
	for( int i=0; i < TempCollection.size() ;i++ )
	{
		for( int j=0; j < TempCollection[i].signature.size() ;j++ )
			if ( TempCollection[i].signature[j] != antigen_signature[j] ) break;
			else if( j == TempCollection[i].signature.size()-1 ){ TempCollection[i].probability_weight++; return true; }
	}
	return false;
}

void AntigenLibrary::add_antigen(const std::vector<double> &antigen_signature, const double current_time)
{
	if ( check_library(antigen_signature) ) return; // This signature already exists
	else
	{
		antigen AntigenTemp;
		AntigenTemp.signature = antigen_signature;
		AntigenTemp.ID = Collection.size() + TempCollection.size();
		AntigenTemp.time = current_time;
		AntigenTemp.probability_weight = 1;
		TempCollection.push_back(AntigenTemp);
	}
	return;
}

void AntigenLibrary::update_collection(const double current_time, const double intervalToRecord)
{
	std::vector<int> IndexesToRemove;
	for( int i=0; i < TempCollection.size() ;i++ )
	{
		if ( current_time - TempCollection[i].time > intervalToRecord || TempCollection.size() == 1 ) // If have a certain time or it's the first antigen
		{
			Collection.push_back(TempCollection[i]);
			IndexesToRemove.push_back(i);
		}
	}
	for( int i=0; i < IndexesToRemove.size() ;i++ )
	{
		TempCollection[ IndexesToRemove[i] ] = TempCollection[TempCollection.size()-1]; // last element replace this element on vector
		TempCollection.pop_back(); // remove the last element from vecotr
	}


	return;
}

std::vector<double> AntigenLibrary::sample_antigen()
{
	std::vector<int> prob(Collection.size(),0);
	for(int i=0; i < Collection.size(); i++) prob[i] = Collection[i].probability_weight;
	std::discrete_distribution<int> discrete_dist(prob.begin(),prob.end());
	int index = discrete_dist(generator);
	return Collection[index].signature;
}

void AntigenLibrary::print()
{
	std::cout << "\n----------------- Neoantigen Collection ----------------" << std::endl;
	for( int i=0; i < Collection.size() ;i++ )
	{
		std::cout << "Antigen ID: " << Collection[i].ID << " Time: " << Collection[i].time << " Weight: " << Collection[i].probability_weight << " Signature: [ ";
		for( int j=0; j < Collection[i].signature.size() ;j++ )
			std::cout << Collection[i].signature[j] << " ";
		std::cout << "]" << std::endl;
	}

	std::cout << "-------------- Waiting Neoantigen Collection ----------------" << std::endl;
	for( int i=0; i < TempCollection.size() ;i++ )
	{
		std::cout << "Antigen ID: " << TempCollection[i].ID << " Time: " << TempCollection[i].time << " Weight: " << TempCollection[i].probability_weight << " Signature: [ ";
		for( int j=0; j < TempCollection[i].signature.size() ;j++ )
			std::cout << TempCollection[i].signature[j] << " ";
		std::cout << "]\n" << std::endl;
	}
}

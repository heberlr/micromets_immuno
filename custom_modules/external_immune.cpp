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
	extern double TH1;
	extern double TH2;
	extern double TCt;
	extern double Tht;
	extern double GridCOUNT;
	static double dC = parameters.doubles( "TC_death_rate" );
	static double pT1 = parameters.doubles( "max_activation_TC" );
	static double pT2 = parameters.doubles( "half_max_activation_TC" );
	static double dT1 = parameters.doubles( "max_clearance_TC" );
	static double dT2 = parameters.doubles( "half_max_clearance_TC" );
	static double Tc0 = parameters.doubles( "TC_population_threshold" );
	static double immunevolume = 1;
	static double dDm = parameters.doubles( "DM_decay" );
	static double sTh1 = 7e-4;
	static double pTh1 = 8.3e-6; //new value 1.5e-5;
	static double dTh1 = 7e-7;
	static double mTh = 1.5e-5; // new value 1.56e-5
	static double sTh2 = 3e-5 ; // new value 2.8e-5
	static double pTh2 = 2e-6;
	static double ro = 1;
	static double CD8_Tcell_recruitment_rate = parameters.doubles( "T_Cell_Recruitment" );

	//double lypmh_scale = GridCOUNT / 3e6; // default 5e6
	double lypmh_scale = 0.02; // old version

	// actual model goes here

	double x[4][6]={{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};//initialize x
	double f[4][6]={{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};//initialize f
	int j;

	// TC update
	double dR_TC = dC * Tc0;

	/* 	// DM Tc recruitment
	double dR_TCD = pT1 * C[0]/immunevolume * C[1]/immunevolume / ( C[0]/immunevolume + pT2/immunevolume) ;

	// DM Tc decay
	double dR_TC16 = dT1 * C[0]/immunevolume * C[1]/immunevolume / ( C[0]/immunevolume + dT2/immunevolume) ;

	// TC decay
	double dR_TC14 = dC * C[1] / immunevolume ;

	// DM decay
	double dR_DM = dDm * C[0] / immunevolume; */




	x[0][0] = DM/lypmh_scale;
	x[0][1] = TC; //initial values
	x[0][2] = TH1; //initial values
	x[0][3] = TH2; //initial values
	x[0][4] = TCt/lypmh_scale;
	x[0][5] = Tht/lypmh_scale;

	for(j = 0; j < 4; j++){
		f[j][0] = {-dDm*x[j][0]}; //define function
		f[j][1] = {dR_TC-dC*x[j][1]+pT1*x[j][0]*x[j][1]/(x[j][0]+pT2)-dT1*x[j][0]*x[j][1]/(x[j][0]+dT2)};// /* TEST */ + (1e6 - x[j][1])/1e6};
		f[j][2] = {(sTh1*x[j][2])/((1+x[j][3])*(1+x[j][3]))+(pTh1*x[j][0]*x[j][2]*x[j][2])/((1+x[j][3])*(1+x[j][3]))-(dTh1*x[j][0]*x[j][2]*x[j][2]*x[j][2])/(500+x[j][3])-mTh*x[j][2]}; //define function
		f[j][3] = {(sTh2*x[j][3])/(1+x[j][3])+(pTh2*(ro+x[j][2])*x[j][0]*x[j][3]*x[j][3])/((1+x[j][3])*(1+x[j][2]+x[j][3]))-mTh*x[j][3]}; //define function
		f[j][4] = {CD8_Tcell_recruitment_rate*x[j][1]}; //define function
		f[j][5] = {CD8_Tcell_recruitment_rate*(x[j][2]+x[j][3])}; //define function
		if (j== 0 || j==1){
			x[j+1][0]=x[0][0]+dt/2*f[j][0]; //first and second x approximations
			x[j+1][1]=x[0][1]+dt/2*f[j][1]; //first and second x approximations
			x[j+1][2]=x[0][2]+dt/2*f[j][2]; //first and second x approximations
			x[j+1][3]=x[0][3]+dt/2*f[j][3]; //first and second x approximations
			x[j+1][4]=x[0][4]+dt/2*f[j][4]; //first and second x approximations
			x[j+1][5]=x[0][5]+dt/2*f[j][5]; //first and second x approximations
		}
		if (j== 2){
			x[j+1][0]=x[0][0]+dt*f[j][0]; //third approximation
			x[j+1][1]=x[0][1]+dt*f[j][1]; //third approximation
			x[j+1][2]=x[0][2]+dt*f[j][2]; //third approximation
			x[j+1][3]=x[0][3]+dt*f[j][3]; //third approximation
			x[j+1][4]=x[0][4]+dt*f[j][4]; //third approximation
			x[j+1][5]=x[0][5]+dt*f[j][5]; //third approximation
		}
	}

	//std::cout << dt*(f[0][0]/6+f[1][0]/3+f[2][0]/3+f[3][0]/6) << std::endl;

	DM=(x[0][0]+dt*(f[0][0]/6+f[1][0]/3+f[2][0]/3+f[3][0]/6))*lypmh_scale;
	TC=x[0][1]+dt*(f[0][1]/6+f[1][1]/3+f[2][1]/3+f[3][1]/6);
	TH1=x[0][2]+dt*(f[0][2]/6+f[1][2]/3+f[2][2]/3+f[3][2]/6); //detirmine n+1
	TH2=x[0][3]+dt*(f[0][3]/6+f[1][3]/3+f[2][3]/3+f[3][3]/6); //detirmine n+1
	TCt=(x[0][4]+dt*(f[0][4]/6+f[1][4]/3+f[2][4]/3+f[3][4]/6))*lypmh_scale;
	Tht=(x[0][5]+dt*(f[0][5]/6+f[1][5]/3+f[2][5]/3+f[3][5]/6))*lypmh_scale;

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
			else if( j == Collection[i].signature.size()-1 ) return true;
	}
	// Check temporary collection
	for( int i=0; i < TempCollection.size() ;i++ )
	{
		for( int j=0; j < TempCollection[i].signature.size() ;j++ )
			if ( TempCollection[i].signature[j] != antigen_signature[j] ) break;
			else if( j == TempCollection[i].signature.size()-1 ) return true;
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
		AntigenTemp.ID = Collection.size();
		AntigenTemp.time = current_time;
		TempCollection.push_back(AntigenTemp);
	}
	return;
}

void RemoveAntigenFromTempCollection( const int index)
{

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
	std::uniform_int_distribution<int> unif_discrete_dist(0,Collection.size()-1);
	int index = unif_discrete_dist(generator);
	return Collection[index].signature;
}

void AntigenLibrary::print()
{
	std::cout << "\n----------------- Neoantigen Collection ----------------" << std::endl;
	for( int i=0; i < Collection.size() ;i++ )
	{
		std::cout << "Antigen ID: " << Collection[i].ID << " Time: " << Collection[i].time << " Signature: [ ";
		for( int j=0; j < Collection[i].signature.size() ;j++ )
			std::cout << Collection[i].signature[j] << " ";
		std::cout << "]" << std::endl;
	}

	std::cout << "-------------- Waiting Neoantigen Collection ----------------" << std::endl;
	for( int i=0; i < TempCollection.size() ;i++ )
	{
		std::cout << "Antigen ID: " << TempCollection[i].ID << " Time: " << TempCollection[i].time << " Signature: [ ";
		for( int j=0; j < TempCollection[i].signature.size() ;j++ )
			std::cout << TempCollection[i].signature[j] << " ";
		std::cout << "]\n" << std::endl;
	}
}

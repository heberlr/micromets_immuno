#include "./external_immune.h"

using namespace PhysiCell;

std::string external_immune_version = "0.1.0";

Submodel_Information external_immune_info;

// generator sample from distribution
extern std::default_random_engine generator;

extern LymphNode lympNode;

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

	external_immune_info.register_model();

	return;
}

void external_immune_model( double dt )
{
	static double dC = parameters.doubles( "TC_death_rate" );
	static double pT1 = parameters.doubles( "activation_rate_TC" );
	static double pT2 = parameters.doubles( "half_max_activation_TC" );
	static double dT1 = parameters.doubles( "clearance_rate_TC" );
	static double dT2 = parameters.doubles( "half_max_clearance_TC" );
	static double Tc0 = parameters.doubles( "TC_population_threshold" );
	static double dDm = parameters.doubles( "DM_decay" );
	static double sTh1 = parameters.doubles( "Th1_max_activation" );
	static double pTh1 = parameters.doubles( "Th1_damping" ); //8.3e-6;
	static double dTh1 = parameters.doubles( "Th1_decay" );
	static double mTh = parameters.doubles( "Th_base_decay" ); //1.5e-5;
	static double sTh2 = parameters.doubles( "Th2_self_feeback" ); //3e-5
	static double pTh2 = parameters.doubles( "Th2_max_conversion" );
	static double ro = parameters.doubles( "Th1_Th2_conversion_weight" );
	static double CD8_Tcell_recruitment_rate = parameters.doubles( "T_Cell_Recruitment" );

	double lypmh_scale = lympNode.GridCOUNT / 5e5;

	// actual model goes here
	double x[4][6]={{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};//initialize x
	double f[4][6]={{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};//initialize f
	int j;

	// TC update
	double dR_TC = dC * Tc0;

	extern std::vector<int>history;

	x[0][0] = (lympNode.DM+history.back())/lypmh_scale;
	x[0][1] = lympNode.TC; //initial values
	x[0][2] = lympNode.TH1; //initial values
	x[0][3] = lympNode.TH2; //initial values
	x[0][4] = lympNode.TCt/lypmh_scale;
	x[0][5] = lympNode.Tht/lypmh_scale;

	for(j = 0; j < 4; j++){
		f[j][0] = {-dDm*x[j][0]}; //define function
		f[j][1] = {dR_TC-dC*x[j][1]+pT1*((1000000-x[j][1])/(1000000))*x[j][0]*x[j][1]/(x[j][0]+pT2)-dT1*x[j][0]*x[j][1]/(x[j][0]+dT2)};
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

	lympNode.DM=(x[0][0]+dt*(f[0][0]/6+f[1][0]/3+f[2][0]/3+f[3][0]/6))*lypmh_scale;
	lympNode.TC=x[0][1]+dt*(f[0][1]/6+f[1][1]/3+f[2][1]/3+f[3][1]/6);
	lympNode.TH1=x[0][2]+dt*(f[0][2]/6+f[1][2]/3+f[2][2]/3+f[3][2]/6); //detirmine n+1
	lympNode.TH2=x[0][3]+dt*(f[0][3]/6+f[1][3]/3+f[2][3]/3+f[3][3]/6); //detirmine n+1
	lympNode.TCt=(x[0][4]+dt*(f[0][4]/6+f[1][4]/3+f[2][4]/3+f[3][4]/6))*lypmh_scale;
	lympNode.Tht=(x[0][5]+dt*(f[0][5]/6+f[1][5]/3+f[2][5]/3+f[3][5]/6))*lypmh_scale;


	return;
}

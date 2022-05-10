#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

#include "./submodel_data_structures.h"

#ifndef __external_immune__
#define __external_immune__

extern Submodel_Information external_immune_info;

void external_immune_model_setup( void );
void external_immune_model( double dt );

class LymphNode{
	public:
		double DCAMOUNT = 0.0;
		double GridCOUNT = 1.0;
		double DM;
		double TC;
		double TH1;
		double TH2;
		double TCt;
		double Tht;
		void InitialCondition(double dm,double tc,double th1,double th2,double tct,double tht){
			DM = dm; TC=tc; TH1=th1; TH2=th2; TCt=tct; Tht=tht;
		}
};

#endif

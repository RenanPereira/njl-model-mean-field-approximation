#ifndef NJLDIMENSIONFULCOUPLINGS_H
#define NJLDIMENSIONFULCOUPLINGS_H

#include <vector>
#include "generalPhysicsAndMath.h"

using namespace std;


enum lagrangianInteractions { interactions_4SP_det,
							  interactions_4SP_det_4VP,
							  interactions_4SP_det_4VP_4VIPI,
							  interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI,
							  interactions_4SP_det_4VP_8VP_8SPVP,
							  interactions_4SP_det_8SP,
							  interactions_4SP_det_8SP_4VP_8VP,
							  interactions_4SP_det_8SP_4VP_8VP_8SPVP,
							  interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI,
							  interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI };


class NJLDimensionfulCouplings
{
//quark interactions: SP=ScalarPseudoscalar , VP=VectorPseudovector , VIPI=VectorIsovectorPseudovectorIsovector
private:
	lagrangianInteractions interactions;

	//4 quark interaction couplings[GeV-2]
	double fourQuarkSPCoupling = 0.0;//g
	double fourQuarkVPCoupling = 0.0;//gOmega
	double fourQuarkVIPICoupling = 0.0;//gRho
	
	//tHooft determinant quark interaction couplings[GeV-5 for SU3]
	double determinantCoupling = 0.0;//kappa
	
	//8 quark interaction couplings[GeV-8]
	double eightQuarkSPOziViolatingCoupling = 0.0;//g1
	double eightQuarkSPNonOziViolatingCoupling = 0.0;//g2
	double eightQuarkVPCoupling = 0.0;//gOmegaOmega
	double eightQuarkVIPICoupling = 0.0;//gRhoRho
	double eightQuarkVPVIPICoupling = 0.0;//gOmegaRho
	double eightQuarkSPVPCoupling = 0.0;//gSigmaOmega
	double eightQuarkSPVIPICoupling = 0.0;//gSigmaRho

public:
	NJLDimensionfulCouplings(){};
	NJLDimensionfulCouplings(lagrangianInteractions , double , double );
	NJLDimensionfulCouplings(lagrangianInteractions , double , double , double );
	NJLDimensionfulCouplings(lagrangianInteractions , double , double , double , double );
	NJLDimensionfulCouplings(lagrangianInteractions , double , double , double , double , double );
	NJLDimensionfulCouplings(lagrangianInteractions , double , double , double , double , double , double );
	NJLDimensionfulCouplings(lagrangianInteractions , double , double , double , double , double , double , double );
	NJLDimensionfulCouplings(lagrangianInteractions , double , double , double , double , double , double , double , double , double );
	NJLDimensionfulCouplings(lagrangianInteractions , double , double , double , double , double , double , double , double , double , double , double );

	lagrangianInteractions getLagrangianInteractions(){ return interactions; };
	double getFourQuarkSPCoupling(){ return fourQuarkSPCoupling; };
	double getFourQuarkVPCoupling(){ return fourQuarkVPCoupling; };
	double getFourQuarkVIPICoupling(){ return fourQuarkVIPICoupling; };
	double getDeterminantCoupling(){ return determinantCoupling; };
	double getEightQuarkSPOziViolatingCoupling(){ return eightQuarkSPOziViolatingCoupling; };
	double getEightQuarkSPNonOziViolatingCoupling(){ return eightQuarkSPNonOziViolatingCoupling; };
	double getEightQuarkVPCoupling(){ return eightQuarkVPCoupling; };
	double getEightQuarkVIPICoupling(){ return eightQuarkVIPICoupling; };
	double getEightQuarkVPVIPICoupling(){ return eightQuarkVPVIPICoupling; };
	double getEightQuarkSPVPCoupling(){ return eightQuarkSPVPCoupling; };
	double getEightQuarkSPVIPICoupling(){ return eightQuarkSPVIPICoupling; };
};



#endif
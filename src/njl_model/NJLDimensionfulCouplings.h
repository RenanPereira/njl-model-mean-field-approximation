#ifndef NJLDIMENSIONFULCOUPLINGS_H
#define NJLDIMENSIONFULCOUPLINGS_H

#include <vector>
#include <string>
#include <map>
#include "ini_file_parser/IniFileParser.h"

using namespace std;


enum LagrangianInteractions { 
	interactions_4SP_det,
	interactions_4SP_det_4VP,
	interactions_4SP_det_4VP_8VP,
	interactions_4SP_det_4VP_8VP_12VP,
	interactions_4SP_det_4VP_8VP_12VP_16VP,
	interactions_4SP_det_4VP_4VIPI,
	interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI,
	interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI,
	interactions_4SP_det_4VP_8VP_8SPVP,
	interactions_4SP_det_8SP,
	interactions_4SP_det_8SP_4VP_8VP,
	interactions_4SP_det_8SP_4VP_8VP_8SPVP,
	interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI,
	interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI,
	interactions_4SP_det_multiVP
};

// Create the mapping between the enum and its string representation
static const map<LagrangianInteractions, string> LagrangianInteractionsMap = 
{
    {interactions_4SP_det, "interactions_4SP_det"},
    {interactions_4SP_det_4VP, "interactions_4SP_det_4VP"},
    {interactions_4SP_det_4VP_8VP, "interactions_4SP_det_4VP_8VP"},
    {interactions_4SP_det_4VP_8VP_12VP, "interactions_4SP_det_4VP_8VP_12VP"},
    {interactions_4SP_det_4VP_8VP_12VP_16VP, "interactions_4SP_det_4VP_8VP_12VP_16VP"},
    {interactions_4SP_det_4VP_4VIPI, "interactions_4SP_det_4VP_4VIPI"},
    {interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI, "interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI"},
    {interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI, "interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI"},
    {interactions_4SP_det_4VP_8VP_8SPVP, "interactions_4SP_det_4VP_8VP_8SPVP"},
    {interactions_4SP_det_8SP, "interactions_4SP_det_8SP"},
    {interactions_4SP_det_8SP_4VP_8VP, "interactions_4SP_det_8SP_4VP_8VP"},
    {interactions_4SP_det_8SP_4VP_8VP_8SPVP, "interactions_4SP_det_8SP_4VP_8VP_8SPVP"},
    {interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI, "interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI"},
    {interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI, "interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI"},
    {interactions_4SP_det_multiVP, "interactions_4SP_det_multiVP"}
};

string toStringLagrangianInteractions(LagrangianInteractions );

LagrangianInteractions stringToLagrangianInteractions(const string& );

bool isValidLagrangianInteractions(const std::string& );


class NJLDimensionfulCouplings
{
//quark interactions: SP=ScalarPseudoscalar , VP=VectorPseudovector , VIPI=VectorIsovectorPseudovectorIsovector
private:
	LagrangianInteractions interactions;

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

	//12 quark interaction couplings[GeV-14]
	double twelveQuarkVPCoupling = 0.0;//gOmega3

	//16 quark interaction couplings[GeV-20]
	double sixteenQuarkVPCoupling = 0.0;//gOmega4

	//multi VP quark interaction couplings: the number of elements corresponds to the number of considered increasing VP interactions
	vector<double> multiQuarkVPCoupling = {};
	bool interactionsIncludeMultiQuarkVPCouplings = false;

public:
	NJLDimensionfulCouplings(){};
	NJLDimensionfulCouplings(LagrangianInteractions , double , double );
	NJLDimensionfulCouplings(LagrangianInteractions , double , double , double );
	NJLDimensionfulCouplings(LagrangianInteractions , double , double , double , double );
	NJLDimensionfulCouplings(LagrangianInteractions , double , double , double , double , double );
	NJLDimensionfulCouplings(LagrangianInteractions , double , double , double , double , double , double );
	NJLDimensionfulCouplings(LagrangianInteractions , double , double , double , double , double , double , double );
	NJLDimensionfulCouplings(LagrangianInteractions , double , double , double , double , double , double , double , double , double );
	NJLDimensionfulCouplings(LagrangianInteractions , double , double , double , double , double , double , double , double , double , double , double );
	NJLDimensionfulCouplings(LagrangianInteractions , double , double , vector<double> );

	LagrangianInteractions getLagrangianInteractions(){ return interactions; };
	
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

	double getTwelveQuarkVPCoupling(){ return twelveQuarkVPCoupling; };

	double getSixteenQuarkVPCoupling(){ return sixteenQuarkVPCoupling; };

	vector<double> getMultiQuarkVPCoupling(){ return multiQuarkVPCoupling; }
	double getMultiQuarkVPCoupling(int i){ return multiQuarkVPCoupling[i]; }
	bool getInteractionsIncludeMultiQuarkVPCouplings(){ return interactionsIncludeMultiQuarkVPCouplings; }
	int numberOfMultiQuarkVPCoupling(){ return int( multiQuarkVPCoupling.size() ); }

	void errorWrongConstructor();
};


vector<double> multiQuarkVPCouplingWithDimensions(vector<double> , double );

bool validateNJLDimensionfulCouplings(const IniFileParser& , string , string );

#endif
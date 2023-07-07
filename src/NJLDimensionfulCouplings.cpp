#include <cmath>
#include <iostream>
#include "NJLDimensionfulCouplings.h"

using namespace std;


void NJLDimensionfulCouplings::errorWrongConstructor()
{
	cout << "Wrong constructor for this lagrangian interaction!\n"; abort();
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2)
{	
	interactions = interactionsAux;

	if ( interactions==interactions_4SP_det )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3)
{
	interactions = interactionsAux;

	if ( interactions==interactions_4SP_det_4VP )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4)
{	
	interactions = interactionsAux;

	if ( interactions==interactions_4SP_det_8SP )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
	}
	else if ( interactions==interactions_4SP_det_4VP_4VIPI )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		fourQuarkVIPICoupling = c4;
	}
	else if ( interactions==interactions_4SP_det_4VP_8VP )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		eightQuarkVPCoupling = c4;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5)
{
	interactions = interactionsAux;

	if ( interactions==interactions_4SP_det_4VP_8VP_8SPVP )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		eightQuarkVPCoupling = c4;
		eightQuarkSPVPCoupling = c5;
	}
	else if ( interactions==interactions_4SP_det_4VP_8VP_12VP )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		eightQuarkVPCoupling = c4;
		twelveQuarkVPCoupling = c5;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6)
{
	interactions = interactionsAux;

	if ( interactions==interactions_4SP_det_8SP_4VP_8VP )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
		fourQuarkVPCoupling = c5;
		eightQuarkVPCoupling = c6;
	}
	else if ( interactions==interactions_4SP_det_4VP_8VP_12VP_16VP )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		eightQuarkVPCoupling = c4;
		twelveQuarkVPCoupling = c5;
		sixteenQuarkVPCoupling = c6;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6, double c7)
{
	interactions = interactionsAux;

	if ( interactions==interactions_4SP_det_8SP_4VP_8VP_8SPVP )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
		fourQuarkVPCoupling = c5;
		eightQuarkVPCoupling = c6;
		eightQuarkSPVPCoupling = c7;
	}
	else if( interactions==interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		fourQuarkVIPICoupling = c4;
		eightQuarkVPCoupling = c5;
		eightQuarkVIPICoupling = c6;
		eightQuarkVPVIPICoupling = c7;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8, double c9)
{	
	interactions = interactionsAux;

	if ( interactions==interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
		fourQuarkVPCoupling = c5;
		fourQuarkVIPICoupling = c6;
		eightQuarkVPCoupling = c7;
		eightQuarkVIPICoupling = c8;
		eightQuarkVPVIPICoupling = c9;
	}
	else if ( interactions==interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		fourQuarkVIPICoupling = c4;
		eightQuarkVPCoupling = c5;
		eightQuarkVIPICoupling = c6;
		eightQuarkVPVIPICoupling = c7;
		eightQuarkSPVPCoupling = c8;
		eightQuarkSPVIPICoupling = c9;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8, double c9, double c10, double c11)
{
	interactions = interactionsAux;

	if ( interactions==interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
		fourQuarkVPCoupling = c5;
		fourQuarkVIPICoupling = c6;
		eightQuarkVPCoupling = c7;
		eightQuarkVIPICoupling = c8;
		eightQuarkVPVIPICoupling = c9;
		eightQuarkSPVPCoupling = c10;
		eightQuarkSPVIPICoupling = c11;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, vector<double> v)
{	
	interactions = interactionsAux;

	if ( interactions==interactions_4SP_det_multiVP )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		multiQuarkVPCoupling = v;
		if ( int( multiQuarkVPCoupling.size() )!=0 )
		{
			interactionsIncludeMultiQuarkVPCouplings = true;
		}
	}
	else{ errorWrongConstructor(); }
}


vector<double> multiQuarkVPCouplingWithDimensions(vector<double> multiQuarkVPCouplingWithoutDimensions, double couplingGeVMinus2)
{
	vector<double> multiQuarkVPDimensionfullCouplings;

    for (int i = 0; i < int( multiQuarkVPCouplingWithoutDimensions.size() ); ++i)
    {	
    	int p = 1 + 3*i;
    	double gs = couplingGeVMinus2;
    	double ci = multiQuarkVPCouplingWithoutDimensions[i];
    	multiQuarkVPDimensionfullCouplings.push_back( ci*pow(gs,p) );
    }

    return multiQuarkVPDimensionfullCouplings;
}


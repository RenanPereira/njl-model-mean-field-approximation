#include <cmath>
#include <iostream>
#include "NJLDimensionfulCouplings.h"

using namespace std;


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2)
{	
	if ( interactionsAux==interactions_4SP_det )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
	}
	else{ cout << "Wrong constructor for this lagrangian interaction!\n"; abort(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3)
{
	if ( interactionsAux==interactions_4SP_det_4VP )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
	}
	else{ cout << "Wrong constructor for this lagrangian interaction!\n"; abort(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4)
{
	if ( interactionsAux==interactions_4SP_det_8SP )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
	}
	else if ( interactionsAux==interactions_4SP_det_4VP_4VIPI )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		fourQuarkVIPICoupling = c4;
	}
	else{ cout << "Wrong constructor for this lagrangian interaction!\n"; abort(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5)
{
	if ( interactionsAux==interactions_4SP_det_4VP_8VP_8SPVP )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		eightQuarkVPCoupling = c4;
		eightQuarkSPVPCoupling = c5;
	}
	else{ cout << "Wrong constructor for this lagrangian interaction!\n"; abort(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6)
{
	if ( interactionsAux==interactions_4SP_det_8SP_4VP_8VP )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
		fourQuarkVPCoupling = c5;
		eightQuarkVPCoupling = c6;
	}
	else{ cout << "Wrong constructor for this lagrangian interaction!\n"; abort(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6, double c7)
{
	if ( interactionsAux==interactions_4SP_det_8SP_4VP_8VP_8SPVP )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
		fourQuarkVPCoupling = c5;
		eightQuarkVPCoupling = c6;
		eightQuarkSPVPCoupling = c7;
	}
	else if( interactionsAux==interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		fourQuarkVIPICoupling = c4;
		eightQuarkVPCoupling = c5;
		eightQuarkVIPICoupling = c6;
		eightQuarkVPVIPICoupling = c7;
	}
	else{ cout << "Wrong constructor for this lagrangian interaction!\n"; abort(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8, double c9)
{
	if ( interactionsAux==interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI )
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
	else{ cout << "Wrong constructor for this lagrangian interaction!\n"; abort(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(lagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8, double c9, double c10, double c11)
{
	if ( interactionsAux==interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI )
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
	else{ cout << "Wrong constructor for this lagrangian interaction!\n"; abort(); }
}


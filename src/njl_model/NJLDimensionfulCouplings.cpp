#include <cmath>
#include <iostream>
#include "njl_model/NJLDimensionfulCouplings.h"

using namespace std;


string toStringLagrangianInteractions(lagrangianInteractions interaction) 
{
    switch (interaction) 
	{
        case interactions_4SP_det:
            return "interactions_4SP_det";
        case interactions_4SP_det_4VP:
            return "interactions_4SP_det_4VP";
        case interactions_4SP_det_4VP_8VP:
            return "interactions_4SP_det_4VP_8VP";
        case interactions_4SP_det_4VP_8VP_12VP:
            return "interactions_4SP_det_4VP_8VP_12VP";
        case interactions_4SP_det_4VP_8VP_12VP_16VP:
            return "interactions_4SP_det_4VP_8VP_12VP_16VP";
        case interactions_4SP_det_4VP_4VIPI:
            return "interactions_4SP_det_4VP_4VIPI";
        case interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI:
            return "interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI";
        case interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI:
            return "interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI";
        case interactions_4SP_det_4VP_8VP_8SPVP:
            return "interactions_4SP_det_4VP_8VP_8SPVP";
        case interactions_4SP_det_8SP:
            return "interactions_4SP_det_8SP";
        case interactions_4SP_det_8SP_4VP_8VP:
            return "interactions_4SP_det_8SP_4VP_8VP";
        case interactions_4SP_det_8SP_4VP_8VP_8SPVP:
            return "interactions_4SP_det_8SP_4VP_8VP_8SPVP";
        case interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI:
            return "interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI";
        case interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI:
            return "interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI";
        case interactions_4SP_det_multiVP:
            return "interactions_4SP_det_multiVP";
        default:
			cout << "Invalid lagrangianInteractions value, returning Unknown" << endl;
            return "Unknown";
    }
}


lagrangianInteractions fromStringLagrangianInteractions(const string& interactionStr) 
{
    if (interactionStr == "interactions_4SP_det") 
	{
        return interactions_4SP_det;
    } 
	else if (interactionStr == "interactions_4SP_det_4VP") 
	{
        return interactions_4SP_det_4VP;
    } 
	else if (interactionStr == "interactions_4SP_det_4VP_8VP") 
	{
        return interactions_4SP_det_4VP_8VP;
    } 
	else if (interactionStr == "interactions_4SP_det_4VP_8VP_12VP") 
	{
        return interactions_4SP_det_4VP_8VP_12VP;
    } 
	else if (interactionStr == "interactions_4SP_det_4VP_8VP_12VP_16VP") 
	{
        return interactions_4SP_det_4VP_8VP_12VP_16VP;
    } 
	else if (interactionStr == "interactions_4SP_det_4VP_4VIPI") 
	{
        return interactions_4SP_det_4VP_4VIPI;
    } 
	else if (interactionStr == "interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI") 
	{
        return interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI;
    } 
	else if (interactionStr == "interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI") 
	{
        return interactions_4SP_det_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI;
    } 
	else if (interactionStr == "interactions_4SP_det_4VP_8VP_8SPVP") 
	{
        return interactions_4SP_det_4VP_8VP_8SPVP;
    } 
	else if (interactionStr == "interactions_4SP_det_8SP") 
	{
        return interactions_4SP_det_8SP;
    } 
	else if (interactionStr == "interactions_4SP_det_8SP_4VP_8VP") 
	{
        return interactions_4SP_det_8SP_4VP_8VP;
    } 
	else if (interactionStr == "interactions_4SP_det_8SP_4VP_8VP_8SPVP") 
	{
        return interactions_4SP_det_8SP_4VP_8VP_8SPVP;
    } 
	else if (interactionStr == "interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI") 
	{
        return interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI;
    } 
	else if (interactionStr == "interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI") 
	{
        return interactions_4SP_det_8SP_4VP_4VIPI_8VP_8VIPI_8VPVIPI_8SPVP_8SPVIPI;
    } 
	else if (interactionStr == "interactions_4SP_det_multiVP") 
	{
        return interactions_4SP_det_multiVP;
    } 
	else 
	{
        cout << "Invalid lagrangianInteractions string: " + interactionStr + ". Aborting!\n";
        abort();
    }
}


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


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


NJLDimensionfulCouplings extractSU3NJL3DCutoffDimensionfulCouplings(const IniFileParser& config)
{	
	double cutoff_GeV = config.getDouble("SU3NJL3DCutoffModelParameters", "cutoff_GeV");

    lagrangianInteractions interactionTerms = fromStringLagrangianInteractions(config.getValue("NJLDimensionfulCouplings", "lagrangianInteractions"));
	cout << "interactionTerms = " << toStringLagrangianInteractions(interactionTerms) << endl;

	if ( interactionTerms==interactions_4SP_det )
	{
		double fourQuarkSPCouplingCutoff2 = config.getDouble("NJLDimensionfulCouplings", "fourQuarkSPCouplingCutoff2");
		double determinantCouplingCutoff5 = config.getDouble("NJLDimensionfulCouplings", "determinantCouplingCutoff5");
		
		cout << "fourQuarkSPCouplingCutoff2 = " << fourQuarkSPCouplingCutoff2 << endl;
		cout << "determinantCouplingCutoff5 = " << determinantCouplingCutoff5 << endl;
		
		double gs = fourQuarkSPCouplingCutoff2/pow(cutoff_GeV, 2);
		double kappa = determinantCouplingCutoff5/pow(cutoff_GeV, 5);
		
		NJLDimensionfulCouplings couplingsSU3NJL3DCutoff(interactionTerms, gs, kappa);

		return couplingsSU3NJL3DCutoff;
	}
	else if( interactionTerms==interactions_4SP_det_8SP )
	{
		double fourQuarkSPCouplingCutoff2 = config.getDouble("NJLDimensionfulCouplings", "fourQuarkSPCouplingCutoff2");
		double determinantCouplingCutoff5 = config.getDouble("NJLDimensionfulCouplings", "determinantCouplingCutoff5");
		double eightQuarkSPOziViolatingCouplingCutoff8 = config.getDouble("NJLDimensionfulCouplings", "eightQuarkSPOziViolatingCouplingCutoff8");
		double eightQuarkSPNonOziViolatingCouplingCutoff8 = config.getDouble("NJLDimensionfulCouplings", "eightQuarkSPNonOziViolatingCouplingCutoff8");
		
		cout << "fourQuarkSPCouplingCutoff2 = " << fourQuarkSPCouplingCutoff2 << endl;
		cout << "determinantCouplingCutoff5 = " << determinantCouplingCutoff5 << endl;
		cout << "eightQuarkSPOziViolatingCouplingCutoff8 = " << eightQuarkSPOziViolatingCouplingCutoff8 << endl;
		cout << "eightQuarkSPNonOziViolatingCouplingCutoff8 = " << eightQuarkSPNonOziViolatingCouplingCutoff8 << endl;
		
		double gs = fourQuarkSPCouplingCutoff2/pow(cutoff_GeV, 2);
		double kappa = determinantCouplingCutoff5/pow(cutoff_GeV, 5);
		double g1 = eightQuarkSPOziViolatingCouplingCutoff8/pow(cutoff_GeV, 8);
		double g2 = eightQuarkSPNonOziViolatingCouplingCutoff8/pow(cutoff_GeV, 8);    
		
		NJLDimensionfulCouplings couplingsSU3NJL3DCutoff(interactionTerms, gs, kappa, g1, g2);

		return couplingsSU3NJL3DCutoff;
	}
	else
	{
		cout << "Setting vector lagrangian interactions via ini file feeding is not yet supported. Aborting...";
		abort();
	}
}

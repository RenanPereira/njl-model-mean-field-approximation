#include <cmath>
#include <iostream>
#include "njl_model/NJLDimensionfulCouplings.h"
#include "njl_model/NJLDimensionlessCouplings.h"

using namespace std;


string toString(LagrangianInteractions interaction) 
{
	// Check if the method exists in the map using count
    if (LagrangianInteractionsMap.count(interaction))
    {
        return LagrangianInteractionsMap.at(interaction);
    } 
    else 
    {
        cout << "Error: LagrangianInteractions not found in map! Returning UNKNOWN." << endl;
        return "UNKNOWN";
    }
}


LagrangianInteractions stringToLagrangianInteractions(const string& interactionString) 
{
	// Iterate over the map with explicit type
    for (map<LagrangianInteractions, string>::const_iterator it = LagrangianInteractionsMap.begin(); it != LagrangianInteractionsMap.end(); ++it) 
    {
        if (it->second == interactionString) 
        {
            return it->first;
        }
    }

    cout << "Invalid LagrangianInteractions string: " + interactionString + ". Aborting!\n";
    abort();
}


bool isValidLagrangianInteractions(const string& interactionString)
{
	bool isLagrangianInteractionsValid = false;
    // Iterate over the map with explicit type
    for (map<LagrangianInteractions, string>::const_iterator it = LagrangianInteractionsMap.begin(); it != LagrangianInteractionsMap.end(); ++it) 
    {
        if (it->second == interactionString) 
        {
            isLagrangianInteractionsValid = true;
            break;
        }
    }

    if( isLagrangianInteractionsValid==false )
    {
        cout << "The value " + interactionString + " is not a LagrangianInteractions!\n";
    }

    return isLagrangianInteractionsValid;
}


void NJLDimensionfulCouplings::errorWrongConstructor()
{
	cout << "Wrong constructor for this lagrangian interaction!\n"; abort();
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(LagrangianInteractions interactionsAux, double c1, double c2)
{	
	interactions = interactionsAux;

	if ( interactions==SP4Q_DET2NFQ )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(LagrangianInteractions interactionsAux, double c1, double c2, double c3)
{
	interactions = interactionsAux;

	if ( interactions==SP4Q_DET2NFQ_VP4Q )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(LagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4)
{	
	interactions = interactionsAux;

	if ( interactions==SP4Q_DET2NFQ_SP8Q )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
	}
	else if ( interactions==SP4Q_DET2NFQ_VP4Q_VIPI4Q )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		fourQuarkVIPICoupling = c4;
	}
	else if ( interactions==SP4Q_DET2NFQ_VP4Q_VP8Q )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		eightQuarkVPCoupling = c4;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(LagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5)
{
	interactions = interactionsAux;

	if ( interactions==SP4Q_DET2NFQ_VP4Q_VP8Q_SPVP8Q )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		eightQuarkVPCoupling = c4;
		eightQuarkSPVPCoupling = c5;
	}
	else if ( interactions==SP4Q_DET2NFQ_VP4Q_VP8Q_VP12Q )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		fourQuarkVPCoupling = c3;
		eightQuarkVPCoupling = c4;
		twelveQuarkVPCoupling = c5;
	}
	else{ errorWrongConstructor(); }
}


NJLDimensionfulCouplings::NJLDimensionfulCouplings(LagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6)
{
	interactions = interactionsAux;

	if ( interactions==SP4Q_DET2NFQ_SP8Q_VP4Q_VP8Q )
	{	
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
		fourQuarkVPCoupling = c5;
		eightQuarkVPCoupling = c6;
	}
	else if ( interactions==SP4Q_DET2NFQ_VP4Q_VP8Q_VP12Q_VP16Q )
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


NJLDimensionfulCouplings::NJLDimensionfulCouplings(LagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6, double c7)
{
	interactions = interactionsAux;

	if ( interactions==SP4Q_DET2NFQ_SP8Q_VP4Q_VP8Q_SPVP8Q )
	{
		fourQuarkSPCoupling = c1;
		determinantCoupling = c2;
		eightQuarkSPOziViolatingCoupling = c3;
		eightQuarkSPNonOziViolatingCoupling = c4;
		fourQuarkVPCoupling = c5;
		eightQuarkVPCoupling = c6;
		eightQuarkSPVPCoupling = c7;
	}
	else if( interactions==SP4Q_DET2NFQ_VP4Q_VIPI4Q_VP8Q_VIPI8Q_VPVIPI8Q )
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


NJLDimensionfulCouplings::NJLDimensionfulCouplings(LagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8, double c9)
{	
	interactions = interactionsAux;

	if ( interactions==SP4Q_DET2NFQ_SP8Q_VP4Q_VIPI4Q_VP8Q_VIPI8Q_VPVIPI8Q )
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
	else if ( interactions==SP4Q_DET2NFQ_VP4Q_VIPI4Q_VP8Q_VIPI8Q_VPVIPI8Q_SPVP8Q_SPVIPI8Q )
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


NJLDimensionfulCouplings::NJLDimensionfulCouplings(LagrangianInteractions interactionsAux, double c1, double c2, double c3, double c4, double c5, double c6, double c7, double c8, double c9, double c10, double c11)
{
	interactions = interactionsAux;

	if ( interactions==SP4Q_DET2NFQ_SP8Q_VP4Q_VIPI4Q_VP8Q_VIPI8Q_VPVIPI8Q_SPVP8Q_SPVIPI8Q )
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


NJLDimensionfulCouplings::NJLDimensionfulCouplings(LagrangianInteractions interactionsAux, double c1, double c2, vector<double> v)
{	
	interactions = interactionsAux;

	if ( interactions==SP4Q_DET2NFQ_VPMULTIQ )
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


bool validateNJLDimensionfulCouplings(const IniFileParser& config, string sectionNJLDimensionfulCouplings, string keyLagrangianInteractions)
{    
    LagrangianInteractions interaction = stringToLagrangianInteractions(config.getValue(sectionNJLDimensionfulCouplings, keyLagrangianInteractions));
	
    if ( interaction==SP4Q_DET2NFQ )
	{   
        bool fourQuarkSPCouplingPresent = config.isKeyPresent(sectionNJLDimensionfulCouplings, NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING);
        bool determinantCouplingPresent = config.isKeyPresent(sectionNJLDimensionfulCouplings, NJLDimensionlessCouplings::DETERMINANT_COUPLING);

        if( fourQuarkSPCouplingPresent && 
            determinantCouplingPresent )
        {
            return true;
        }
        else{ return false; }
	}
	else if( interaction==SP4Q_DET2NFQ_SP8Q )
	{   
        bool fourQuarkSPCouplingPresent = config.isKeyPresent(sectionNJLDimensionfulCouplings, NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING);
        bool determinantCouplingPresent = config.isKeyPresent(sectionNJLDimensionfulCouplings, NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING);
        bool eightQuarkSPOziViolatingCouplingPresent = config.isKeyPresent(sectionNJLDimensionfulCouplings, NJLDimensionlessCouplings::EIGHT_QUARK_SP_OZI_COUPLING);
        bool eightQuarkSPNonOziViolatingCouplingPresent = config.isKeyPresent(sectionNJLDimensionfulCouplings, NJLDimensionlessCouplings::EIGHT_QUARK_SP_NON_OZI_COUPLING);

        if( fourQuarkSPCouplingPresent && 
            determinantCouplingPresent &&
            eightQuarkSPOziViolatingCouplingPresent &&
            eightQuarkSPNonOziViolatingCouplingPresent )
        {
            return true;
        }
        else{ return false; }
	}
	else
	{
		cout << "Setting vector lagrangian interactions via ini file feeding is not yet supported.";
        return false;
	}
}

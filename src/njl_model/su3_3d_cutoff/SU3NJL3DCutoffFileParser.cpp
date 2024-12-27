#include <iostream>
#include <cmath>
#include "gsl_wrapper/root_solver_gsl.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "njl_model/njl_regularization_schemes.h"
#include "njl_model/NJLDimensionlessCouplings.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCalculator.h"


using namespace std;


NJLDimensionfulCouplings SU3NJL3DCutoffVacuumFileParser::extractDimensionfulCouplings(const IniFileParser& config)
{	
	double cutoffInGeV = config.getDouble(SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
                                          SU3NJL3DCutoffConfigKeys::ModelParameters::cutoffInGeV);

    LagrangianInteractions interaction = 
        stringToLagrangianInteractions(
            config.getValue(SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
                            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::lagrangianInteractions)
        );
	cout << "interaction = " << toStringLagrangianInteractions(interaction) << endl;

	if ( interaction==interactions_4SP_det )
	{
		double fourQuarkSPCouplingCutoff2 = 
            config.getDouble(SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
                             NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING);
		double determinantCouplingCutoff5 = 
            config.getDouble(SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
                             NJLDimensionlessCouplings::DETERMINANT_COUPLING);
		
		cout << NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING + " = " << fourQuarkSPCouplingCutoff2 << endl;
		cout << NJLDimensionlessCouplings::DETERMINANT_COUPLING + " = " << determinantCouplingCutoff5 << endl;
		
		double gs = fourQuarkSPCouplingCutoff2/pow(cutoffInGeV, 2);
		double kappa = determinantCouplingCutoff5/pow(cutoffInGeV, 5);
		
		NJLDimensionfulCouplings couplingsSU3NJL3DCutoff(interaction, gs, kappa);

		return couplingsSU3NJL3DCutoff;
	}
	else if( interaction==interactions_4SP_det_8SP )
	{
		double fourQuarkSPCouplingCutoff2 = 
            config.getDouble(SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
                             NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING);
		double determinantCouplingCutoff5 = 
            config.getDouble(SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
                             NJLDimensionlessCouplings::DETERMINANT_COUPLING);
		double eightQuarkSPOziViolatingCouplingCutoff8 = 
            config.getDouble(SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
                             NJLDimensionlessCouplings::EIGHT_QUARK_SP_OZI_COUPLING);
		double eightQuarkSPNonOziViolatingCouplingCutoff8 = 
            config.getDouble(SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
                             NJLDimensionlessCouplings::EIGHT_QUARK_SP_NON_OZI_COUPLING);
		
		cout << NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING << " = " << fourQuarkSPCouplingCutoff2 << endl;
		cout << NJLDimensionlessCouplings::DETERMINANT_COUPLING + " = " << determinantCouplingCutoff5 << endl;
		cout << NJLDimensionlessCouplings::EIGHT_QUARK_SP_OZI_COUPLING + " = " << eightQuarkSPOziViolatingCouplingCutoff8 << endl;
		cout << NJLDimensionlessCouplings::EIGHT_QUARK_SP_NON_OZI_COUPLING + " = " << eightQuarkSPNonOziViolatingCouplingCutoff8 << endl;
		
		double gs = fourQuarkSPCouplingCutoff2/pow(cutoffInGeV, 2);
		double kappa = determinantCouplingCutoff5/pow(cutoffInGeV, 5);
		double g1 = eightQuarkSPOziViolatingCouplingCutoff8/pow(cutoffInGeV, 8);
		double g2 = eightQuarkSPNonOziViolatingCouplingCutoff8/pow(cutoffInGeV, 8);    
		
		NJLDimensionfulCouplings couplingsSU3NJL3DCutoff(interaction, gs, kappa, g1, g2);

		return couplingsSU3NJL3DCutoff;
	}
	else
	{
		cout << "Setting vector lagrangian interactions via ini file feeding is not yet supported. Aborting...";
		abort();
	}
}


bool SU3NJL3DCutoffVacuumFileParser::validateFileQuality() const 
{
    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = {
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
        SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::section
    };  
    for (int i = 0; i < int(requiredSections.size()); ++i) 
    {
        string section = requiredSections[i];
        if (config.getSectionsData(section).empty()) 
        {   
            allRequiredSectionsPresent = false;
            cout << "Missing required section: " << section << endl;
        }
    }

    // Validate individual sections
    bool areModelParametersValid = validateModelParameters();
    bool areDimensionfulCouplingsValid = validateDimensionfulCouplings();
    bool areGapEquationsVacuumParametersValid = validateGapEquationsVacuumParameters();

    // The function return true only if all tests passed
    if ( allRequiredSectionsPresent &&
         areModelParametersValid &&
         areDimensionfulCouplingsValid &&
         areGapEquationsVacuumParametersValid )
    { 
        return true; 
    }
    else{ return false; }
    
}

bool SU3NJL3DCutoffVacuumFileParser::validateModelParameters() const 
{   
    // Validate regularizationScheme
    string regularizationScheme = 
        config.getValue(SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
                        SU3NJL3DCutoffConfigKeys::ModelParameters::regularizationScheme);
    bool isRegularizationSchemeValid = isValidNJL3DCutoffRegularizationScheme(regularizationScheme);
    if( !isRegularizationSchemeValid ){ cout << invalidFileMessage << endl; }

    // Ensure cutoffInGeV>0
    bool isCutoffInGeVValid = 
        config.validatePositiveDouble(SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
                                      SU3NJL3DCutoffConfigKeys::ModelParameters::cutoffInGeV, 
                                      invalidFileMessage, 
                                      SU3NJL3DCutoffConfigKeys::ModelParameters::cutoffInGeV + " > 0 must be satisfied.");

    // Ensure upQuarkCurrentMassInGeV>0
    bool isUpQuarkCurrentMassInGeVValid = 
        config.validatePositiveDouble(SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
                                      SU3NJL3DCutoffConfigKeys::ModelParameters::upQuarkCurrentMassInGeV, 
                                      invalidFileMessage, 
                                      SU3NJL3DCutoffConfigKeys::ModelParameters::upQuarkCurrentMassInGeV + " > 0 must be satisfied.");

    // Ensure downQuarkCurrentMassInGeV>0
    bool isDownQuarkCurrentMassInGeVValid = 
        config.validatePositiveDouble(SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
                                      SU3NJL3DCutoffConfigKeys::ModelParameters::downQuarkCurrentMassInGeV, 
                                      invalidFileMessage, 
                                      SU3NJL3DCutoffConfigKeys::ModelParameters::downQuarkCurrentMassInGeV + " > 0 must be satisfied.");

    // Ensure strangeQuarkCurrentMassInGeV>0
    bool isStrangeQuarkCurrentMassInGeVValid = 
        config.validatePositiveDouble(SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
                                      SU3NJL3DCutoffConfigKeys::ModelParameters::strangeQuarkCurrentMassInGeV, 
                                      invalidFileMessage, 
                                      SU3NJL3DCutoffConfigKeys::ModelParameters::strangeQuarkCurrentMassInGeV + " > 0 must be satisfied.");

    // The function return true only if all tests passed
    if ( isRegularizationSchemeValid && 
         isCutoffInGeVValid &&
         isUpQuarkCurrentMassInGeVValid &&
         isDownQuarkCurrentMassInGeVValid &&
         isStrangeQuarkCurrentMassInGeVValid )
    { 
        return true; 
    }
    else{ return false; }
}

bool SU3NJL3DCutoffVacuumFileParser::validateDimensionfulCouplings() const 
{
    // Validate lagrangianInteractions
    string lagrangianInteractions = config.getValue(SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
                                                    SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::lagrangianInteractions);
    bool isLagrangianInteractionsValid = isValidLagrangianInteractions(lagrangianInteractions);
    if( !isLagrangianInteractionsValid ){ cout << invalidFileMessage << endl; }

    bool areCouplingsValid = false;
    if ( isLagrangianInteractionsValid )
    {
        areCouplingsValid = validateNJLDimensionfulCouplings(config,
                                                             SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,
                                                             SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::lagrangianInteractions);
    }

    return areCouplingsValid;
}

bool SU3NJL3DCutoffVacuumFileParser::validateGapEquationsVacuumParameters() const 
{   
    // Ensure gapPrecision>0
    bool isGapPrecisionValid = 
        config.validatePositiveDouble(SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::section,
                                      SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::gapPrecision, 
                                      invalidFileMessage, 
                                      SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::gapPrecision + " > 0 must be satisfied.");

    // Ensure rootFindingMethod is valid
    string rootFindingMethod = 
        config.getValue(SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::section, 
                        SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::rootFindingMethod);
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(rootFindingMethod);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    bool isUpQuarkMassGuessValid = 
        config.validatePositiveDouble(SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::section, 
                                      SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::upQuarkMassGuess, 
                                      invalidFileMessage, 
                                      SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::upQuarkMassGuess + " > 0 must be satisfied.");

    bool isDownQuarkMassGuessValid = 
        config.validatePositiveDouble(SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::section, 
                                      SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::downQuarkMassGuess, 
                                      invalidFileMessage, 
                                      SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::downQuarkMassGuess + " > 0 must be satisfied.");

    bool isStrangeQuarkMassGuessValid = 
        config.validatePositiveDouble(SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::section, 
                                      SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::strangeQuarkMassGuess, 
                                      invalidFileMessage, 
                                      SU3NJL3DCutoffConfigKeys::GapEquationsVacuumParameters::strangeQuarkMassGuess + " > 0 must be satisfied.");

    // The function return true only if all tests passed
    if ( isGapPrecisionValid && 
         isRootFindingMethodValid &&
         isDownQuarkMassGuessValid &&
         isStrangeQuarkMassGuessValid &&
         isUpQuarkMassGuessValid )
    { 
        return true; 
    }
    else{ return false; }
}


void SU3NJL3DCutoffVacuumFileParser::evaluateVacuumMasses() const
{
    // Get SU3 NJL 3D Cutoff Model Parameters
    cout << "\nSU3NJL3DCutoffModelParameters:" << endl;

    string parameterSetName = config.getValue("SU3NJL3DCutoffModelParameters", "parameterSetName");
    string regularizationScheme = config.getValue("SU3NJL3DCutoffModelParameters", "regularizationScheme");
    double cutoffInGeV = config.getDouble("SU3NJL3DCutoffModelParameters", "cutoffInGeV");
    double upQuarkCurrentMassInGeV = config.getDouble("SU3NJL3DCutoffModelParameters", "upQuarkCurrentMassInGeV");
    double downQuarkCurrentMassInGeV = config.getDouble("SU3NJL3DCutoffModelParameters", "downQuarkCurrentMassInGeV");
    double strangeQuarkCurrentMassInGeV = config.getDouble("SU3NJL3DCutoffModelParameters", "strangeQuarkCurrentMassInGeV");
    
    cout << "parameterSetName = " << parameterSetName << endl;
    cout << "regularizationScheme = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << "cutoffInGeV = " << cutoffInGeV << endl;
    cout << "upQuarkCurrentMassInGeV = " << upQuarkCurrentMassInGeV << endl;
    cout << "downQuarkCurrentMassInGeV = " << downQuarkCurrentMassInGeV << endl;
    cout << "strangeQuarkCurrentMassInGeV = " << strangeQuarkCurrentMassInGeV << endl;


    // Get SU3 NJL 3D Cutoff Dimensionful Couplings
    cout << "\nSU3NJL3DCutoffGapEquationsVacuumParameters:" << endl;
    NJLDimensionfulCouplings couplings = SU3NJL3DCutoffVacuumFileParser::extractDimensionfulCouplings(config);


    // SU3 NJL 3D Cutoff Gap Equations Vacuum Parameters
    cout << "\nSU3NJL3DCutoffGapEquationsVacuumParameters: " << endl;

    double gapPrecision = config.getDouble("SU3NJL3DCutoffGapEquationsVacuumParameters", "gapPrecision");
    string rootFindingMethod = config.getValue("SU3NJL3DCutoffGapEquationsVacuumParameters", "rootFindingMethod");
    double upQuarkMassGuess = config.getDouble("SU3NJL3DCutoffGapEquationsVacuumParameters", "upQuarkMassGuess");
    double downQuarkMassGuess = config.getDouble("SU3NJL3DCutoffGapEquationsVacuumParameters", "downQuarkMassGuess");
    double strangeQuarkMassGuess = config.getDouble("SU3NJL3DCutoffGapEquationsVacuumParameters", "strangeQuarkMassGuess");

    cout << "gapPrecision = " << gapPrecision << endl;
    cout << "rootFindingMethod = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(rootFindingMethod)) << endl;
    cout << "upQuarkMassGuess = " << upQuarkMassGuess << endl;
    cout << "downQuarkMassGuess = " << downQuarkMassGuess << endl;
    cout << "strangeQuarkMassGuess = " << strangeQuarkMassGuess << endl;


    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
                                        cutoffInGeV, 
                                        couplings, 
                                        upQuarkCurrentMassInGeV, 
                                        downQuarkCurrentMassInGeV, 
                                        strangeQuarkCurrentMassInGeV);
    parameters.setParameterSetName(parameterSetName);

    evaluateSU3NJL3DCutoffVacuumMasses(parameters, 
                                       gapPrecision, 
                                       stringToMultiRootFindingMethod(rootFindingMethod), 
                                       upQuarkMassGuess, 
                                       downQuarkMassGuess, 
                                       strangeQuarkMassGuess);
}

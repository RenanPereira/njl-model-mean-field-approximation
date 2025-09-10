#include <iostream>
#include <cmath>
#include "gsl_wrapper/root_solver_gsl.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "njl_model/njl_regularization_schemes.h"
#include "njl_model/NJLDimensionlessCouplings.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedTempRhoBEqualChemPot.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.h"

using namespace std;


NJLDimensionfulCouplings SU3NJL3DCutoffVacuumFileParser::extractDimensionfulCouplings(const IniFileParser& config)
{	
	double cutoff = config.getDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section,                                   
        SU3NJL3DCutoffConfigKeys::ModelParameters::cutoff
    );

    LagrangianInteractions interaction = stringToLagrangianInteractions(
        config.getValue(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,             
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::lagrangianInteractions
        )
    );
	cout << "interaction = " << toStringLagrangianInteractions(interaction) << endl;

	if ( interaction==SP4Q_DET2NFQ )
	{
		double fourQuarkSPCouplingCutoff2 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,              
            NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING
        );
		double determinantCouplingCutoff5 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,              
            NJLDimensionlessCouplings::DETERMINANT_COUPLING
        );
		
		cout << NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING + " = " << fourQuarkSPCouplingCutoff2 << endl;
		cout << NJLDimensionlessCouplings::DETERMINANT_COUPLING + " = " << determinantCouplingCutoff5 << endl;
		
		double gs = fourQuarkSPCouplingCutoff2/pow(cutoff, 2);
		double kappa = determinantCouplingCutoff5/pow(cutoff, 5);
		
		NJLDimensionfulCouplings couplingsSU3NJL3DCutoff(interaction, gs, kappa);

		return couplingsSU3NJL3DCutoff;
	}
	else if( interaction==SP4Q_DET2NFQ_SP8Q )
	{
		double fourQuarkSPCouplingCutoff2 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING
        );
		double determinantCouplingCutoff5 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::DETERMINANT_COUPLING
        );
		double eightQuarkSPOziViolatingCouplingCutoff8 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::EIGHT_QUARK_SP_OZI_COUPLING
        );
		double eightQuarkSPNonOziViolatingCouplingCutoff8 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::EIGHT_QUARK_SP_NON_OZI_COUPLING
        );
		
		cout << NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING << " = " << fourQuarkSPCouplingCutoff2 << endl;
		cout << NJLDimensionlessCouplings::DETERMINANT_COUPLING + " = " << determinantCouplingCutoff5 << endl;
		cout << NJLDimensionlessCouplings::EIGHT_QUARK_SP_OZI_COUPLING + " = " << eightQuarkSPOziViolatingCouplingCutoff8 << endl;
		cout << NJLDimensionlessCouplings::EIGHT_QUARK_SP_NON_OZI_COUPLING + " = " << eightQuarkSPNonOziViolatingCouplingCutoff8 << endl;
		
		double gs = fourQuarkSPCouplingCutoff2/pow(cutoff, 2);
		double kappa = determinantCouplingCutoff5/pow(cutoff, 5);
		double g1 = eightQuarkSPOziViolatingCouplingCutoff8/pow(cutoff, 8);
		double g2 = eightQuarkSPNonOziViolatingCouplingCutoff8/pow(cutoff, 8);    
		
		NJLDimensionfulCouplings couplingsSU3NJL3DCutoff(interaction, gs, kappa, g1, g2);

		return couplingsSU3NJL3DCutoff;
	}
	else
	{
		cout << "Setting vector lagrangian interactions via ini file feeding is not yet supported. Aborting...";
		abort();
	}
}


bool SU3NJL3DCutoffVacuumFileParser::validateFileQualityEvaluateVacuumMasses() const 
{
    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section
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
    bool areVacuumMassesParametersValid = validateVacuumMassesParameters();

    // The function return true only if all tests passed
    if (allRequiredSectionsPresent &&
        areModelParametersValid &&
        areDimensionfulCouplingsValid &&
        areVacuumMassesParametersValid)
    { 
        return true; 
    }
    else{ return false; }
}

bool SU3NJL3DCutoffVacuumFileParser::validateModelParameters() const 
{   
    // Validate regularizationScheme
    string regularizationScheme = config.getValue(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::regularizationScheme
    );
    bool isRegularizationSchemeValid = isValidNJL3DCutoffRegularizationScheme(regularizationScheme);
    if( !isRegularizationSchemeValid ){ cout << invalidFileMessage << endl; }

    // Ensure cutoff>0
    bool isCutoffValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::cutoff, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::cutoff + " > 0 must be satisfied."
    );

    // Ensure upQuarkCurrentMass>0
    bool isUpQuarkCurrentMassValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::upQuarkCurrentMass, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::upQuarkCurrentMass + " > 0 must be satisfied."
    );

    // Ensure downQuarkCurrentMassInGeV>0
    bool isDownQuarkCurrentMassValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::downQuarkCurrentMass, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::downQuarkCurrentMass + " > 0 must be satisfied."
    );

    // Ensure strangeQuarkCurrentMass>0
    bool isStrangeQuarkCurrentMassValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::strangeQuarkCurrentMass, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::strangeQuarkCurrentMass + " > 0 must be satisfied."
    );

    // The function return true only if all tests passed
    if (isRegularizationSchemeValid && 
        isCutoffValid &&
        isUpQuarkCurrentMassValid &&
        isDownQuarkCurrentMassValid &&
        isStrangeQuarkCurrentMassValid)
    { 
        return true; 
    }
    else{ return false; }
}

bool SU3NJL3DCutoffVacuumFileParser::validateDimensionfulCouplings() const 
{
    // Validate lagrangianInteractions
    string lagrangianInteractions = config.getValue(
        SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,                                             
        SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::lagrangianInteractions
    );
    bool isLagrangianInteractionsValid = isValidLagrangianInteractions(lagrangianInteractions);
    if( !isLagrangianInteractionsValid ){ cout << invalidFileMessage << endl; }

    bool areCouplingsValid = false;
    if ( isLagrangianInteractionsValid )
    {
        areCouplingsValid = validateNJLDimensionfulCouplings(
            config,
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::lagrangianInteractions
        );
    }

    return areCouplingsValid;
}

bool SU3NJL3DCutoffVacuumFileParser::validateVacuumMassesParameters() const 
{   
    // Ensure precisionVacuum>0
    bool isPrecisionVacuumValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section,
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::precisionVacuum, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::precisionVacuum + " > 0 must be satisfied."
    );

    // Ensure methodVacuum is valid
    string methodVacuum = config.getValue(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::methodVacuum
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodVacuum);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // Ensure up quark mass is valid
    bool isUpQuarkMassGuessValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::upQuarkMassGuess, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::upQuarkMassGuess + " > 0 must be satisfied."
    );

    // Ensure down quark mass is valid
    bool isDownQuarkMassGuessValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::downQuarkMassGuess, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::downQuarkMassGuess + " > 0 must be satisfied."
    );

    // Ensure strange quark mass is valid
    bool isStrangeQuarkMassGuessValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::strangeQuarkMassGuess, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::strangeQuarkMassGuess + " > 0 must be satisfied."
    );

    // The function return true only if all tests passed
    if (isPrecisionVacuumValid && 
        isRootFindingMethodValid &&
        isDownQuarkMassGuessValid &&
        isStrangeQuarkMassGuessValid &&
        isUpQuarkMassGuessValid)
    { 
        return true; 
    }
    else{ return false; }
}


void SU3NJL3DCutoffVacuumFileParser::evaluateVacuumMasses() const
{
    // Model Parameters
    namespace MPKeys = SU3NJL3DCutoffConfigKeys::ModelParameters;
    cout << "\n" << MPKeys::section << ": " << endl;

    string parameterSetName = config.getValue(MPKeys::section, MPKeys::parameterSetName);
    string regularizationScheme = config.getValue(MPKeys::section, MPKeys::regularizationScheme);
    double cutoff = config.getDouble(MPKeys::section, MPKeys::cutoff);
    double upQuarkCurrentMass = config.getDouble(MPKeys::section, MPKeys::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MPKeys::section, MPKeys::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MPKeys::section, MPKeys::strangeQuarkCurrentMass);
    
    cout << MPKeys::parameterSetName << " = " << parameterSetName << endl;
    cout << MPKeys::regularizationScheme << " = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MPKeys::cutoff << " = " << cutoff << endl;
    cout << MPKeys::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MPKeys::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MPKeys::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings(config);

    // VacuumMasses
    namespace VMPKeys = SU3NJL3DCutoffConfigKeys::VacuumMassesParameters;
    cout << "\n" << VMPKeys::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMPKeys::section, VMPKeys::precisionVacuum);
    string methodVacuum = config.getValue(VMPKeys::section, VMPKeys::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMPKeys::section, VMPKeys::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMPKeys::section, VMPKeys::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMPKeys::section, VMPKeys::strangeQuarkMassGuess);

    cout << VMPKeys::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMPKeys::methodVacuum << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMPKeys::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMPKeys::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMPKeys::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(
        stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
        cutoff, 
        couplings, 
        upQuarkCurrentMass, 
        downQuarkCurrentMass, 
        strangeQuarkCurrentMass
    );
    parameters.setParameterSetName(parameterSetName);

    SU3NJL3DCutoffVacuum::evaluateVacuumMasses(
        parameters, 
        precisionVacuum, 
        stringToMultiRootFindingMethod(methodVacuum), 
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess
    );
}

bool SU3NJL3DCutoffFixedTempRhoBEqualChemPotFileParser::validateVacuumToFiniteBaryonDensityParameters() const 
{   
    namespace VFBDPKeys = SU3NJL3DCutoffConfigKeys::VacuumToFiniteBaryonDensityParameters;

    // Ensure minimumBaryonDensity>0
    bool isMinimumBaryonDensityValid = config.validatePositiveDouble(
        VFBDPKeys::section,
        VFBDPKeys::minimumBaryonDensity,
        invalidFileMessage + " Invalid value found in section " + VFBDPKeys::section + ".", 
        VFBDPKeys::minimumBaryonDensity + " > 0 must be satisfied."
    );

    // Ensure maximumBaryonDensity>0
    bool isMaximumBaryonDensityValid = config.validatePositiveDouble(
        VFBDPKeys::section,
        VFBDPKeys::maximumBaryonDensity,
        invalidFileMessage + " Invalid value found in section " + VFBDPKeys::section + ".", 
        VFBDPKeys::maximumBaryonDensity + " > 0 must be satisfied."
    );
    
    // Ensure isNumberOfPointsValid>0
    bool isNumberOfPointsValid = config.validatePositiveInteger(
        VFBDPKeys::section,
        VFBDPKeys::numberOfPoints, 
        invalidFileMessage + " Invalid value found in section " + VFBDPKeys::section + ".", 
        VFBDPKeys::numberOfPoints + " > 0 must be satisfied."
    );

    // Ensure precisionZeroTempSol>0
    bool isPrecisionZeroTempSolValid = config.validatePositiveDouble(
        VFBDPKeys::section,
        VFBDPKeys::precisionZeroTempSol, 
        invalidFileMessage + " Invalid value found in section " + VFBDPKeys::section + ".", 
        VFBDPKeys::precisionZeroTempSol + " > 0 must be satisfied."
    );

    // Ensure methodZeroTempSol is valid
    string methodZeroTempSol = config.getValue(
        VFBDPKeys::section, 
        VFBDPKeys::methodZeroTempSol
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodZeroTempSol);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // The function return true only if all tests passed
    if (isMinimumBaryonDensityValid && 
        isMaximumBaryonDensityValid &&
        isNumberOfPointsValid &&
        isPrecisionZeroTempSolValid &&
        isRootFindingMethodValid)
    { 
        return true; 
    }
    else{ return false; }
}

bool SU3NJL3DCutoffFixedTempRhoBEqualChemPotFileParser::validateFirstOrderLineParameters() const 
{   
    namespace FOLPKeys = SU3NJL3DCutoffConfigKeys::FirstOrderLineParameters;

    // Ensure precisionTransitionPointSol>0
    bool isPrecisionTransitionPointSolValid = config.validatePositiveDouble(
        FOLPKeys::section,
        FOLPKeys::precisionTransitionPointSol,
        invalidFileMessage + " Invalid value found in section " + FOLPKeys::section + ".", 
        FOLPKeys::precisionTransitionPointSol + " > 0 must be satisfied."
    );

    // Ensure methodTransitionPointSol is valid
    string methodTransitionPointSol = config.getValue(
        FOLPKeys::section, 
        FOLPKeys::methodTransitionPointSol
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodTransitionPointSol);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // Ensure deltaT>0
    bool isDeltaTValid = config.validatePositiveDouble(
        FOLPKeys::section,
        FOLPKeys::deltaT,
        invalidFileMessage + " Invalid value found in section " + FOLPKeys::section + ".", 
        FOLPKeys::deltaT + " > 0 must be satisfied."
    );

    // Ensure massDifferenceCEP>0
    bool isMassDifferenceCEPValid = config.validatePositiveDouble(
        FOLPKeys::section,
        FOLPKeys::massDifferenceCEP,
        invalidFileMessage + " Invalid value found in section " + FOLPKeys::section + ".", 
        FOLPKeys::massDifferenceCEP + " > 0 must be satisfied."
    );

    // The function return true only if all tests passed
    if (isPrecisionTransitionPointSolValid && 
        isRootFindingMethodValid &&
        isDeltaTValid &&
        isMassDifferenceCEPValid)
    { 
        return true; 
    }
    else{ return false; }
}

bool SU3NJL3DCutoffFixedTempRhoBEqualChemPotFileParser::validateFileQualityEvaluateFirstOrderLine() const
{   
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
    bool vacuumValidations = validateFileQualityEvaluateVacuumMasses();

    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffConfigKeys::VacuumToFiniteBaryonDensityParameters::section,
        SU3NJL3DCutoffConfigKeys::FirstOrderLineParameters::section
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
    bool areVacuumToFiniteBaryonDensityParametersValid = validateVacuumToFiniteBaryonDensityParameters();
    bool areFirstOrderLineParametersValid = validateFirstOrderLineParameters();

    // The function return true only if all tests passed
    if (vacuumValidations && 
        allRequiredSectionsPresent &&
        areVacuumToFiniteBaryonDensityParametersValid &&
        areFirstOrderLineParametersValid)
    { 
        return true; 
    }
    else{ return false; }
}

void SU3NJL3DCutoffFixedTempRhoBEqualChemPotFileParser::evaluateFirstOrderLine() const
{
    // Model Parameters
    namespace MPKeys = SU3NJL3DCutoffConfigKeys::ModelParameters;
    cout << "\n" << MPKeys::section << ": " << endl;

    string parameterSetName = config.getValue(MPKeys::section, MPKeys::parameterSetName);
    string regularizationScheme = config.getValue(MPKeys::section, MPKeys::regularizationScheme);
    double cutoff = config.getDouble(MPKeys::section, MPKeys::cutoff);
    double upQuarkCurrentMass = config.getDouble(MPKeys::section, MPKeys::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MPKeys::section, MPKeys::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MPKeys::section, MPKeys::strangeQuarkCurrentMass);
    
    cout << MPKeys::parameterSetName << " = " << parameterSetName << endl;
    cout << MPKeys::regularizationScheme << " = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MPKeys::cutoff << " = " << cutoff << endl;
    cout << MPKeys::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MPKeys::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MPKeys::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings(config);

    // VacuumMasses
    namespace VMPKeys = SU3NJL3DCutoffConfigKeys::VacuumMassesParameters;
    cout << "\n" << VMPKeys::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMPKeys::section, VMPKeys::precisionVacuum);
    string methodVacuum = config.getValue(VMPKeys::section, VMPKeys::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMPKeys::section, VMPKeys::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMPKeys::section, VMPKeys::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMPKeys::section, VMPKeys::strangeQuarkMassGuess);

    cout << VMPKeys::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMPKeys::methodVacuum << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMPKeys::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMPKeys::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMPKeys::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;

    // VacuumToFiniteBaryonDensity
    namespace VFBDPKeys = SU3NJL3DCutoffConfigKeys::VacuumToFiniteBaryonDensityParameters;
    cout << "\n" << VFBDPKeys::section << ": " << endl;
    
    double minimumBaryonDensity = config.getDouble(VFBDPKeys::section, VFBDPKeys::minimumBaryonDensity);
    double maximumBaryonDensity = config.getDouble(VFBDPKeys::section, VFBDPKeys::maximumBaryonDensity);
    int numberOfPoints = config.getInt(VFBDPKeys::section, VFBDPKeys::numberOfPoints);
    double precisionZeroTempSol = config.getDouble(VFBDPKeys::section, VFBDPKeys::precisionZeroTempSol);
    string methodZeroTempSol = config.getValue(VFBDPKeys::section, VFBDPKeys::methodZeroTempSol);

    cout << VFBDPKeys::minimumBaryonDensity << " = " << minimumBaryonDensity << endl;
    cout << VFBDPKeys::maximumBaryonDensity << " = " << maximumBaryonDensity << endl;
    cout << VFBDPKeys::numberOfPoints << " = " << numberOfPoints << endl;
    cout << VFBDPKeys::precisionZeroTempSol << " = " << precisionZeroTempSol << endl;
    cout << VFBDPKeys::methodZeroTempSol << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodZeroTempSol)) << endl;

    // FirstOrderLine
    namespace FOLPKeys = SU3NJL3DCutoffConfigKeys::FirstOrderLineParameters;
    cout << "\n" << FOLPKeys::section << ":\n";
    
    double precisionTransitionPointSol = config.getDouble(FOLPKeys::section, FOLPKeys::precisionTransitionPointSol);
    string methodTransitionPointSol = config.getValue(FOLPKeys::section, FOLPKeys::methodTransitionPointSol);
    double deltaT = config.getDouble(FOLPKeys::section, FOLPKeys::deltaT);
    double massDifferenceCEP = config.getDouble(FOLPKeys::section, FOLPKeys::massDifferenceCEP);
    
    cout << FOLPKeys::precisionTransitionPointSol  << " = " << precisionTransitionPointSol << endl;
    cout << FOLPKeys::methodTransitionPointSol << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodTransitionPointSol)) << endl;
    cout << FOLPKeys::deltaT << " = " << deltaT << endl;
    cout << FOLPKeys::massDifferenceCEP << " = " << massDifferenceCEP << endl;

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(
        stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
        cutoff, 
        couplings, 
        upQuarkCurrentMass, 
        downQuarkCurrentMass, 
        strangeQuarkCurrentMass
    );
    parameters.setParameterSetName(parameterSetName);

    // Perform calculation
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot::evaluateFirstOrderLine(
        parameters, 
        precisionVacuum, 
        stringToMultiRootFindingMethod(methodVacuum), 
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess, 
        minimumBaryonDensity, 
        maximumBaryonDensity, 
        numberOfPoints, 
        precisionZeroTempSol, 
        stringToMultiRootFindingMethod(methodZeroTempSol), 
        true, 
        precisionTransitionPointSol, 
        stringToMultiRootFindingMethod(methodTransitionPointSol), 
        deltaT, 
        massDifferenceCEP
    );
}

bool SU3NJL3DCutoffFixedChemPotTempCrossSectionsFileParser::validateVacuumToFiniteTemperatureAtZeroChemicalPotentialParameters() const
{
    namespace VFTZCPPKeys = SU3NJL3DCutoffConfigKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;

    // Ensure temperature>0
    bool isTemperatureValid = config.validatePositiveDouble(
        VFTZCPPKeys::section,
        VFTZCPPKeys::temperature,
        invalidFileMessage + " Invalid value found in section " + VFTZCPPKeys::section + ".", 
        VFTZCPPKeys::temperature + " > 0 must be satisfied."
    );

    // Ensure numberOfPointsFromVacToFinTemp>0
    bool isNumberOfPointsFromVacToFinTempValid = config.validatePositiveInteger(
        VFTZCPPKeys::section,
        VFTZCPPKeys::numberOfPointsFromVacToFinTemp, 
        invalidFileMessage + " Invalid value found in section " + VFTZCPPKeys::section + ".", 
        VFTZCPPKeys::numberOfPointsFromVacToFinTemp + " > 0 must be satisfied."
    );

    // Ensure precisionVacToFinTemp>0
    bool isPrecisionVacToFinTempValid = config.validatePositiveDouble(
        VFTZCPPKeys::section,
        VFTZCPPKeys::precisionVacToFinTemp,
        invalidFileMessage + " Invalid value found in section " + VFTZCPPKeys::section + ".", 
        VFTZCPPKeys::precisionVacToFinTemp + " > 0 must be satisfied."
    );

    // Ensure methodVacToFinTemp is valid
    string methodVacToFinTemp = config.getValue(
        VFTZCPPKeys::section, 
        VFTZCPPKeys::methodVacToFinTemp
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodVacToFinTemp);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // The function return true only if all tests passed
    return isTemperatureValid && 
           isNumberOfPointsFromVacToFinTempValid &&
           isPrecisionVacToFinTempValid &&
           isRootFindingMethodValid;
}

bool SU3NJL3DCutoffFixedChemPotTempCrossSectionsFileParser::validateFiniteTemperatureToFiniteChemicalPotentialParameters() const
{
    namespace FTFCPPKeys = SU3NJL3DCutoffConfigKeys::FiniteTemperatureToFiniteChemicalPotentialParameters;

    // Ensure chemPot>=0
    bool isChemPotValid = config.validateNonNegativeDouble(
        FTFCPPKeys::section,
        FTFCPPKeys::chemPot,
        invalidFileMessage + " Invalid value found in section " + FTFCPPKeys::section + ".", 
        FTFCPPKeys::chemPot + " > 0 must be satisfied."
    );

    // Ensure numberOfPointsFromFinTempToFinChemPot>0
    bool isNumberOfPointsFromFinTempToFinChemPotValid = config.validatePositiveInteger(
        FTFCPPKeys::section,
        FTFCPPKeys::numberOfPointsFromFinTempToFinChemPot, 
        invalidFileMessage + " Invalid value found in section " + FTFCPPKeys::section + ".", 
        FTFCPPKeys::numberOfPointsFromFinTempToFinChemPot + " > 0 must be satisfied."
    );

    // Ensure precisionFinTempToFinChemPot>0
    bool isPrecisionFinTempToFinChemPotValid = config.validatePositiveDouble(
        FTFCPPKeys::section,
        FTFCPPKeys::precisionFinTempToFinChemPot,
        invalidFileMessage + " Invalid value found in section " + FTFCPPKeys::section + ".", 
        FTFCPPKeys::precisionFinTempToFinChemPot + " > 0 must be satisfied."
    );

    // Ensure methodFinTempToFinChemPot is valid
    string methodFinTempToFinChemPot = config.getValue(
        FTFCPPKeys::section, 
        FTFCPPKeys::methodFinTempToFinChemPot
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodFinTempToFinChemPot);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // The function return true only if all tests passed
    return isChemPotValid && 
           isNumberOfPointsFromFinTempToFinChemPotValid &&
           isPrecisionFinTempToFinChemPotValid &&
           isRootFindingMethodValid;
}

bool SU3NJL3DCutoffFixedChemPotTempCrossSectionsFileParser::validateCrossSectionsParameters() const
{
    namespace CSPKeys = SU3NJL3DCutoffConfigKeys::CrossSectionsParameters;

    // Ensure propagatorIntegralPrecision>0
    bool isPropagatorIntegralPrecisionValid = config.validatePositiveDouble(
        CSPKeys::section,
        CSPKeys::propagatorIntegralPrecision,
        invalidFileMessage + " Invalid value found in section " + CSPKeys::section + ".", 
        CSPKeys::propagatorIntegralPrecision + " > 0 must be satisfied."
    );

    // Ensure precisionCrossSections>0
    bool isPrecisionCrossSectionsValid = config.validatePositiveDouble(
        CSPKeys::section,
        CSPKeys::precisionCrossSections,
        invalidFileMessage + " Invalid value found in section " + CSPKeys::section + ".", 
        CSPKeys::precisionCrossSections + " > 0 must be satisfied."
    );

    // Ensure numberOfPointsCrossSections>0
    bool isNumberOfPointsCrossSectionsValid = config.validatePositiveInteger(
        CSPKeys::section,
        CSPKeys::numberOfPointsCrossSections, 
        invalidFileMessage + " Invalid value found in section " + CSPKeys::section + ".", 
        CSPKeys::numberOfPointsCrossSections + " > 0 must be satisfied."
    );

    // Ensure numberOfThreads>0
    bool isNumberOfThreadsValid = config.validatePositiveInteger(
        CSPKeys::section,
        CSPKeys::numberOfThreads, 
        invalidFileMessage + " Invalid value found in section " + CSPKeys::section + ".", 
        CSPKeys::numberOfThreads + " > 0 must be satisfied."
    );

    // The function return true only if all tests passed
    return isPropagatorIntegralPrecisionValid && 
           isPrecisionCrossSectionsValid &&
           isNumberOfPointsCrossSectionsValid &&
           isNumberOfThreadsValid;
}

bool SU3NJL3DCutoffFixedChemPotTempCrossSectionsFileParser::validateFileQualityEvaluateCrossSectionsEqualLightMasses() const
{
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
    bool vacuumValidations = validateFileQualityEvaluateVacuumMasses();

    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffConfigKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters::section,
        SU3NJL3DCutoffConfigKeys::FiniteTemperatureToFiniteChemicalPotentialParameters::section,
        SU3NJL3DCutoffConfigKeys::CrossSectionsParameters::section
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
    bool areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid = validateVacuumToFiniteTemperatureAtZeroChemicalPotentialParameters();
    bool areFiniteTemperatureToFiniteChemicalPotentialParametersValid = validateFiniteTemperatureToFiniteChemicalPotentialParameters();
    bool areCrossSectionsParametersValid = validateCrossSectionsParameters();

    // The function return true only if all tests passed
    if (vacuumValidations && 
        allRequiredSectionsPresent && 
        areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid && 
        areFiniteTemperatureToFiniteChemicalPotentialParametersValid &&
        areCrossSectionsParametersValid )
    { 
        return true; 
    }
    else{ return false; }
}

void SU3NJL3DCutoffFixedChemPotTempCrossSectionsFileParser::evaluateCrossSectionsEqualLightMasses() const
{   
    // Model Parameters
    namespace MPKeys = SU3NJL3DCutoffConfigKeys::ModelParameters;
    cout << "\n" << MPKeys::section << ": " << endl;

    string parameterSetName = config.getValue(MPKeys::section, MPKeys::parameterSetName);
    string regularizationScheme = config.getValue(MPKeys::section, MPKeys::regularizationScheme);
    double cutoff = config.getDouble(MPKeys::section, MPKeys::cutoff);
    double upQuarkCurrentMass = config.getDouble(MPKeys::section, MPKeys::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MPKeys::section, MPKeys::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MPKeys::section, MPKeys::strangeQuarkCurrentMass);

    cout << MPKeys::parameterSetName << " = " << parameterSetName << endl;
    cout << MPKeys::regularizationScheme << " = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MPKeys::cutoff << " = " << cutoff << endl;
    cout << MPKeys::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MPKeys::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MPKeys::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings(config);

    // VacuumMassesParameters
    namespace VMPKeys = SU3NJL3DCutoffConfigKeys::VacuumMassesParameters;
    cout << "\n" << VMPKeys::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMPKeys::section, VMPKeys::precisionVacuum);
    string methodVacuum = config.getValue(VMPKeys::section, VMPKeys::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMPKeys::section, VMPKeys::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMPKeys::section, VMPKeys::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMPKeys::section, VMPKeys::strangeQuarkMassGuess);

    cout << VMPKeys::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMPKeys::methodVacuum << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMPKeys::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMPKeys::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMPKeys::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;

    // Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(
        stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
        cutoff, 
        couplings, 
        upQuarkCurrentMass, 
        downQuarkCurrentMass, 
        strangeQuarkCurrentMass
    );
    parameters.setParameterSetName(parameterSetName);

    // VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters
    namespace VFTZCPPKeys = SU3NJL3DCutoffConfigKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    cout << "\n" << VFTZCPPKeys::section << ": " << endl;

    double temperature = config.getDouble(VFTZCPPKeys::section, VFTZCPPKeys::temperature);
    int numberOfPointsFromVacToFinTemp = config.getInt(VFTZCPPKeys::section, VFTZCPPKeys::numberOfPointsFromVacToFinTemp);
    double precisionVacToFinTemp = config.getDouble(VFTZCPPKeys::section, VFTZCPPKeys::precisionVacToFinTemp);
    string methodVacToFinTemp = config.getValue(VFTZCPPKeys::section, VFTZCPPKeys::methodVacToFinTemp);

    cout << VFTZCPPKeys::temperature << " = " << temperature << endl;
    cout << VFTZCPPKeys::numberOfPointsFromVacToFinTemp << " = " << numberOfPointsFromVacToFinTemp << endl;
    cout << VFTZCPPKeys::precisionVacToFinTemp << " = " << precisionVacToFinTemp << endl;
    cout << VFTZCPPKeys::methodVacToFinTemp << " = " << methodVacToFinTemp << endl;

    // FiniteTemperatureToFiniteChemicalPotentialParameters
    namespace FTFCPPKeys = SU3NJL3DCutoffConfigKeys::FiniteTemperatureToFiniteChemicalPotentialParameters;
    cout << "\n" << FTFCPPKeys::section << ": " << endl;
    
    double chemPot = config.getDouble(FTFCPPKeys::section, FTFCPPKeys::chemPot);
    int numberOfPointsFromFinTempToFinChemPot = config.getInt(FTFCPPKeys::section, FTFCPPKeys::numberOfPointsFromFinTempToFinChemPot);
    double precisionFinTempToFinChemPot = config.getDouble(FTFCPPKeys::section, FTFCPPKeys::precisionFinTempToFinChemPot);
    string methodFinTempToFinChemPot = config.getValue(FTFCPPKeys::section, FTFCPPKeys::methodFinTempToFinChemPot);

    cout << FTFCPPKeys::chemPot << " = " << chemPot << endl;
    cout << FTFCPPKeys::numberOfPointsFromFinTempToFinChemPot << " = " << numberOfPointsFromFinTempToFinChemPot << endl;
    cout << FTFCPPKeys::precisionFinTempToFinChemPot << " = " << precisionFinTempToFinChemPot << endl;
    cout << FTFCPPKeys::methodFinTempToFinChemPot << " = " << methodFinTempToFinChemPot << endl;

    // CrossSectionsParameters
    namespace CSPKeys = SU3NJL3DCutoffConfigKeys::CrossSectionsParameters;
    cout << "\n" << CSPKeys::section << ": " << endl;

    double propagatorIntegralPrecision = config.getDouble(CSPKeys::section, CSPKeys::propagatorIntegralPrecision);
    bool largeAngleScatteringContribution = config.getBool(CSPKeys::section, CSPKeys::largeAngleScatteringContribution);
    double precisionCrossSections = config.getDouble(CSPKeys::section, CSPKeys::precisionCrossSections);
    int numberOfPointsCrossSections = config.getInt(CSPKeys::section, CSPKeys::numberOfPointsCrossSections);
    int numberOfThreads = config.getInt(CSPKeys::section, CSPKeys::numberOfThreads);

    cout << CSPKeys::propagatorIntegralPrecision << " = " << propagatorIntegralPrecision << endl;
    cout << CSPKeys::largeAngleScatteringContribution << " = " << largeAngleScatteringContribution << endl;
    cout << CSPKeys::precisionCrossSections << " = " << precisionCrossSections << endl;
    cout << CSPKeys::numberOfPointsCrossSections << " = " << numberOfPointsCrossSections << endl;
    cout << CSPKeys::numberOfThreads << " = " << numberOfThreads << endl;

    // Calculate Cross Sections
    SU3NJL3DCutoffFixedChemPotTemp::evaluateCrossSectionsEqualLightMasses(
    	parameters, 
		precisionVacuum, 
		stringToMultiRootFindingMethod(methodVacuum), 
		upQuarkMassGuess, 
		strangeQuarkMassGuess,
        temperature,
        numberOfPointsFromVacToFinTemp,
        precisionVacToFinTemp,
        stringToMultiRootFindingMethod(methodVacToFinTemp),
        chemPot,
        numberOfPointsFromFinTempToFinChemPot,
        precisionFinTempToFinChemPot,
        stringToMultiRootFindingMethod(methodFinTempToFinChemPot),
        propagatorIntegralPrecision,
        largeAngleScatteringContribution,
        precisionCrossSections,
        numberOfPointsCrossSections,
        numberOfThreads
    );
}


//////////////////////////////////////////////////////////////////////////////////////////////
// On-going refactor: the methods above will be substituted by the ones below

NJLDimensionfulCouplings SU3NJL3DCutoffFileParser::Common::extractDimensionfulCouplings() const
{	
	double cutoff = config.getDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section,                                   
        SU3NJL3DCutoffConfigKeys::ModelParameters::cutoff
    );

    LagrangianInteractions interaction = stringToLagrangianInteractions(
        config.getValue(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,             
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::lagrangianInteractions
        )
    );
	cout << "interaction = " << toStringLagrangianInteractions(interaction) << endl;

	if ( interaction==SP4Q_DET2NFQ )
	{
		double fourQuarkSPCouplingCutoff2 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,              
            NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING
        );
		double determinantCouplingCutoff5 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,              
            NJLDimensionlessCouplings::DETERMINANT_COUPLING
        );
		
		cout << NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING + " = " << fourQuarkSPCouplingCutoff2 << endl;
		cout << NJLDimensionlessCouplings::DETERMINANT_COUPLING + " = " << determinantCouplingCutoff5 << endl;
		
		double gs = fourQuarkSPCouplingCutoff2/pow(cutoff, 2);
		double kappa = determinantCouplingCutoff5/pow(cutoff, 5);
		
		NJLDimensionfulCouplings couplingsSU3NJL3DCutoff(interaction, gs, kappa);

		return couplingsSU3NJL3DCutoff;
	}
	else if( interaction==SP4Q_DET2NFQ_SP8Q )
	{
		double fourQuarkSPCouplingCutoff2 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING
        );
		double determinantCouplingCutoff5 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::DETERMINANT_COUPLING
        );
		double eightQuarkSPOziViolatingCouplingCutoff8 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::EIGHT_QUARK_SP_OZI_COUPLING
        );
		double eightQuarkSPNonOziViolatingCouplingCutoff8 = config.getDouble(
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::EIGHT_QUARK_SP_NON_OZI_COUPLING
        );
		
		cout << NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING << " = " << fourQuarkSPCouplingCutoff2 << endl;
		cout << NJLDimensionlessCouplings::DETERMINANT_COUPLING + " = " << determinantCouplingCutoff5 << endl;
		cout << NJLDimensionlessCouplings::EIGHT_QUARK_SP_OZI_COUPLING + " = " << eightQuarkSPOziViolatingCouplingCutoff8 << endl;
		cout << NJLDimensionlessCouplings::EIGHT_QUARK_SP_NON_OZI_COUPLING + " = " << eightQuarkSPNonOziViolatingCouplingCutoff8 << endl;
		
		double gs = fourQuarkSPCouplingCutoff2/pow(cutoff, 2);
		double kappa = determinantCouplingCutoff5/pow(cutoff, 5);
		double g1 = eightQuarkSPOziViolatingCouplingCutoff8/pow(cutoff, 8);
		double g2 = eightQuarkSPNonOziViolatingCouplingCutoff8/pow(cutoff, 8);    
		
		NJLDimensionfulCouplings couplingsSU3NJL3DCutoff(interaction, gs, kappa, g1, g2);

		return couplingsSU3NJL3DCutoff;
	}
	else
	{
		cout << "Setting vector lagrangian interactions via ini file feeding is not yet supported. Aborting...";
		abort();
	}
}

bool SU3NJL3DCutoffFileParser::Common::validateModelParameters() const
{   
    // Validate regularizationScheme
    string regularizationScheme = config.getValue(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::regularizationScheme
    );
    bool isRegularizationSchemeValid = isValidNJL3DCutoffRegularizationScheme(regularizationScheme);
    if( !isRegularizationSchemeValid ){ cout << invalidFileMessage << endl; }

    // Ensure cutoff>0
    bool isCutoffValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::cutoff, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::cutoff + " > 0 must be satisfied."
    );

    // Ensure upQuarkCurrentMass>0
    bool isUpQuarkCurrentMassValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::upQuarkCurrentMass, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::upQuarkCurrentMass + " > 0 must be satisfied."
    );

    // Ensure downQuarkCurrentMassInGeV>0
    bool isDownQuarkCurrentMassValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::downQuarkCurrentMass, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::downQuarkCurrentMass + " > 0 must be satisfied."
    );

    // Ensure strangeQuarkCurrentMass>0
    bool isStrangeQuarkCurrentMassValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::strangeQuarkCurrentMass, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::ModelParameters::strangeQuarkCurrentMass + " > 0 must be satisfied."
    );

    // The function return true only if all tests passed
    return isRegularizationSchemeValid &&
           isCutoffValid &&
           isUpQuarkCurrentMassValid &&
           isDownQuarkCurrentMassValid &&
           isStrangeQuarkCurrentMassValid;
}

bool SU3NJL3DCutoffFileParser::Common::validateDimensionfulCouplings() const
{
    // Validate lagrangianInteractions
    string lagrangianInteractions = config.getValue(
        SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,                                             
        SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::lagrangianInteractions
    );
    bool isLagrangianInteractionsValid = isValidLagrangianInteractions(lagrangianInteractions);
    if( !isLagrangianInteractionsValid ){ cout << invalidFileMessage << endl; }

    bool areCouplingsValid = false;
    if ( isLagrangianInteractionsValid )
    {
        areCouplingsValid = validateNJLDimensionfulCouplings(
            config,
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section,
            SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::lagrangianInteractions
        );
    }

    return areCouplingsValid;
}

bool SU3NJL3DCutoffFileParser::Common::validateVacuumMassesParameters() const
{   
    // Ensure precisionVacuum>0
    bool isPrecisionVacuumValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section,
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::precisionVacuum, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::precisionVacuum + " > 0 must be satisfied."
    );

    // Ensure methodVacuum is valid
    string methodVacuum = config.getValue(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::methodVacuum
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodVacuum);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // Ensure up quark mass is valid
    bool isUpQuarkMassGuessValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::upQuarkMassGuess, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::upQuarkMassGuess + " > 0 must be satisfied."
    );

    // Ensure down quark mass is valid
    bool isDownQuarkMassGuessValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::downQuarkMassGuess, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::downQuarkMassGuess + " > 0 must be satisfied."
    );

    // Ensure strange quark mass is valid
    bool isStrangeQuarkMassGuessValid = config.validatePositiveDouble(
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::strangeQuarkMassGuess, 
        invalidFileMessage, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::strangeQuarkMassGuess + " > 0 must be satisfied."
    );

    // The function returns true only if all tests passed
    return isPrecisionVacuumValid &&
           isRootFindingMethodValid &&
           isDownQuarkMassGuessValid &&
           isStrangeQuarkMassGuessValid &&
           isUpQuarkMassGuessValid;
}

bool SU3NJL3DCutoffFileParser::Common::validateVacuumToFiniteBaryonDensityParameters() const
{   
    namespace VFBD = SU3NJL3DCutoffConfigKeys::VacuumToFiniteBaryonDensityParameters;

    // Ensure minimumBaryonDensity>0
    bool isMinimumBaryonDensityValid = config.validatePositiveDouble(
        VFBD::section,
        VFBD::minimumBaryonDensity,
        invalidFileMessage + " Invalid value found in section " + VFBD::section + ".", 
        VFBD::minimumBaryonDensity + " > 0 must be satisfied."
    );

    // Ensure maximumBaryonDensity>0
    bool isMaximumBaryonDensityValid = config.validatePositiveDouble(
        VFBD::section,
        VFBD::maximumBaryonDensity,
        invalidFileMessage + " Invalid value found in section " + VFBD::section + ".", 
        VFBD::maximumBaryonDensity + " > 0 must be satisfied."
    );
    
    // Ensure isNumberOfPointsValid>0
    bool isNumberOfPointsValid = config.validatePositiveInteger(
        VFBD::section,
        VFBD::numberOfPoints, 
        invalidFileMessage + " Invalid value found in section " + VFBD::section + ".", 
        VFBD::numberOfPoints + " > 0 must be satisfied."
    );

    // Ensure precisionZeroTempSol>0
    bool isPrecisionZeroTempSolValid = config.validatePositiveDouble(
        VFBD::section,
        VFBD::precisionZeroTempSol, 
        invalidFileMessage + " Invalid value found in section " + VFBD::section + ".", 
        VFBD::precisionZeroTempSol + " > 0 must be satisfied."
    );

    // Ensure methodZeroTempSol is valid
    string methodZeroTempSol = config.getValue(
        VFBD::section, 
        VFBD::methodZeroTempSol
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodZeroTempSol);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // The function return true only if all tests passed
    return isMinimumBaryonDensityValid && 
           isMaximumBaryonDensityValid &&
           isNumberOfPointsValid &&
           isPrecisionZeroTempSolValid &&
           isRootFindingMethodValid;
}

bool SU3NJL3DCutoffFileParser::Common::validateFirstOrderLineParameters() const
{   
    namespace FOLP = SU3NJL3DCutoffConfigKeys::FirstOrderLineParameters;

    // Ensure precisionTransitionPointSol>0
    bool isPrecisionTransitionPointSolValid = config.validatePositiveDouble(
        FOLP::section,
        FOLP::precisionTransitionPointSol,
        invalidFileMessage + " Invalid value found in section " + FOLP::section + ".", 
        FOLP::precisionTransitionPointSol + " > 0 must be satisfied."
    );

    // Ensure methodTransitionPointSol is valid
    string methodTransitionPointSol = config.getValue(
        FOLP::section, 
        FOLP::methodTransitionPointSol
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodTransitionPointSol);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // Ensure deltaT>0
    bool isDeltaTValid = config.validatePositiveDouble(
        FOLP::section,
        FOLP::deltaT,
        invalidFileMessage + " Invalid value found in section " + FOLP::section + ".", 
        FOLP::deltaT + " > 0 must be satisfied."
    );

    // Ensure massDifferenceCEP>0
    bool isMassDifferenceCEPValid = config.validatePositiveDouble(
        FOLP::section,
        FOLP::massDifferenceCEP,
        invalidFileMessage + " Invalid value found in section " + FOLP::section + ".", 
        FOLP::massDifferenceCEP + " > 0 must be satisfied."
    );

    // The function return true only if all tests passed
    return isPrecisionTransitionPointSolValid && 
           isRootFindingMethodValid &&
           isDeltaTValid &&
           isMassDifferenceCEPValid;
}

bool SU3NJL3DCutoffFileParser::Common::validateVacuumToFiniteTemperatureAtZeroChemicalPotentialParameters() const
{
    namespace VFTZCPP = SU3NJL3DCutoffConfigKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;

    // Ensure temperature>0
    bool isTemperatureValid = config.validatePositiveDouble(
        VFTZCPP::section,
        VFTZCPP::temperature,
        invalidFileMessage + " Invalid value found in section " + VFTZCPP::section + ".", 
        VFTZCPP::temperature + " > 0 must be satisfied."
    );

    // Ensure numberOfPointsFromVacToFinTemp>0
    bool isNumberOfPointsFromVacToFinTempValid = config.validatePositiveInteger(
        VFTZCPP::section,
        VFTZCPP::numberOfPointsFromVacToFinTemp, 
        invalidFileMessage + " Invalid value found in section " + VFTZCPP::section + ".", 
        VFTZCPP::numberOfPointsFromVacToFinTemp + " > 0 must be satisfied."
    );

    // Ensure precisionVacToFinTemp>0
    bool isPrecisionVacToFinTempValid = config.validatePositiveDouble(
        VFTZCPP::section,
        VFTZCPP::precisionVacToFinTemp,
        invalidFileMessage + " Invalid value found in section " + VFTZCPP::section + ".", 
        VFTZCPP::precisionVacToFinTemp + " > 0 must be satisfied."
    );

    // Ensure methodVacToFinTemp is valid
    string methodVacToFinTemp = config.getValue(
        VFTZCPP::section, 
        VFTZCPP::methodVacToFinTemp
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodVacToFinTemp);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // The function return true only if all tests passed
    return isTemperatureValid && 
           isNumberOfPointsFromVacToFinTempValid &&
           isPrecisionVacToFinTempValid &&
           isRootFindingMethodValid;
}

bool SU3NJL3DCutoffFileParser::Common::validateLowToHighTemperatureAtZeroChemicalPotentialParameters() const
{   
    namespace LHTZCPP = SU3NJL3DCutoffConfigKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters;

    // Ensure minimumTemp>0
    bool isMinimumTempValid = config.validatePositiveDouble(
        LHTZCPP::section,
        LHTZCPP::minimumTemp,
        invalidFileMessage + " Invalid value found in section " + LHTZCPP::section + ".", 
        LHTZCPP::minimumTemp + " > 0 must be satisfied."
    );

    // Ensure maximumTemp>0
    bool isMaximumTempValid = config.validatePositiveDouble(
        LHTZCPP::section,
        LHTZCPP::maximumTemp,
        invalidFileMessage + " Invalid value found in section " + LHTZCPP::section + ".", 
        LHTZCPP::maximumTemp + " > 0 must be satisfied."
    );

    // Ensure numberOfPointsFromLowToHighTemp>0
    bool isNumberOfPointsFromLowToHighTempValid = config.validatePositiveInteger(
        LHTZCPP::section,
        LHTZCPP::numberOfPointsFromLowToHighTemp, 
        invalidFileMessage + " Invalid value found in section " + LHTZCPP::section + ".", 
        LHTZCPP::numberOfPointsFromLowToHighTemp + " > 0 must be satisfied."
    );

    // Ensure precisionLowToHighTemp>0
    bool isPrecisionLowToHighTempValid = config.validatePositiveDouble(
        LHTZCPP::section,
        LHTZCPP::precisionLowToHighTemp,
        invalidFileMessage + " Invalid value found in section " + LHTZCPP::section + ".", 
        LHTZCPP::precisionLowToHighTemp + " > 0 must be satisfied."
    );

    // Ensure methodLowToHighTemp is valid
    string methodLowToHighTemp = config.getValue(LHTZCPP::section, LHTZCPP::methodLowToHighTemp);
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodLowToHighTemp);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // The function return true only if all tests passed
    return isMinimumTempValid && 
           isMaximumTempValid &&
           isNumberOfPointsFromLowToHighTempValid &&
           isPrecisionLowToHighTempValid &&
           isRootFindingMethodValid;
}

bool SU3NJL3DCutoffFileParser::Common::validateFiniteTemperatureToFiniteChemicalPotentialParameters() const
{
    namespace FTFCPP = SU3NJL3DCutoffConfigKeys::FiniteTemperatureToFiniteChemicalPotentialParameters;

    // Ensure chemPot>=0
    bool isChemPotValid = config.validateNonNegativeDouble(
        FTFCPP::section,
        FTFCPP::chemPot,
        invalidFileMessage + " Invalid value found in section " + FTFCPP::section + ".", 
        FTFCPP::chemPot + " > 0 must be satisfied."
    );

    // Ensure numberOfPointsFromFinTempToFinChemPot>0
    bool isNumberOfPointsFromFinTempToFinChemPotValid = config.validatePositiveInteger(
        FTFCPP::section,
        FTFCPP::numberOfPointsFromFinTempToFinChemPot, 
        invalidFileMessage + " Invalid value found in section " + FTFCPP::section + ".", 
        FTFCPP::numberOfPointsFromFinTempToFinChemPot + " > 0 must be satisfied."
    );

    // Ensure precisionFinTempToFinChemPot>0
    bool isPrecisionFinTempToFinChemPotValid = config.validatePositiveDouble(
        FTFCPP::section,
        FTFCPP::precisionFinTempToFinChemPot,
        invalidFileMessage + " Invalid value found in section " + FTFCPP::section + ".", 
        FTFCPP::precisionFinTempToFinChemPot + " > 0 must be satisfied."
    );

    // Ensure methodFinTempToFinChemPot is valid
    string methodFinTempToFinChemPot = config.getValue(FTFCPP::section, FTFCPP::methodFinTempToFinChemPot);
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodFinTempToFinChemPot);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // The function return true only if all tests passed
    return isChemPotValid && 
           isNumberOfPointsFromFinTempToFinChemPotValid &&
           isPrecisionFinTempToFinChemPotValid &&
           isRootFindingMethodValid;
}

bool SU3NJL3DCutoffFileParser::Common::validateCrossSectionsParameters() const
{
    namespace CSP = SU3NJL3DCutoffConfigKeys::CrossSectionsParameters;

    // Ensure propagatorIntegralPrecision>0
    bool isPropagatorIntegralPrecisionValid = config.validatePositiveDouble(
        CSP::section,
        CSP::propagatorIntegralPrecision,
        invalidFileMessage + " Invalid value found in section " + CSP::section + ".", 
        CSP::propagatorIntegralPrecision + " > 0 must be satisfied."
    );

    // Ensure precisionCrossSections>0
    bool isPrecisionCrossSectionsValid = config.validatePositiveDouble(
        CSP::section,
        CSP::precisionCrossSections,
        invalidFileMessage + " Invalid value found in section " + CSP::section + ".", 
        CSP::precisionCrossSections + " > 0 must be satisfied."
    );

    // Ensure numberOfPointsCrossSections>0
    bool isNumberOfPointsCrossSectionsValid = config.validatePositiveInteger(
        CSP::section,
        CSP::numberOfPointsCrossSections, 
        invalidFileMessage + " Invalid value found in section " + CSP::section + ".", 
        CSP::numberOfPointsCrossSections + " > 0 must be satisfied."
    );

    // Ensure numberOfThreads>0
    bool isNumberOfThreadsValid = config.validatePositiveInteger(
        CSP::section,
        CSP::numberOfThreads, 
        invalidFileMessage + " Invalid value found in section " + CSP::section + ".", 
        CSP::numberOfThreads + " > 0 must be satisfied."
    );

    // The function return true only if all tests passed
    return isPropagatorIntegralPrecisionValid && 
           isPrecisionCrossSectionsValid &&
           isNumberOfPointsCrossSectionsValid &&
           isNumberOfThreadsValid;
}

bool SU3NJL3DCutoffFileParser::Common::validateIntegratedCrossSectionsParameters() const
{
    namespace ICSP = SU3NJL3DCutoffConfigKeys::IntegratedCrossSectionsParameters;

    // Ensure numberOfPointsIntegratedCrossSections>0
    bool isNumberOfPointsIntegratedCrossSectionsValid = config.validatePositiveInteger(
        ICSP::section,
        ICSP::numberOfPointsIntegratedCrossSections, 
        invalidFileMessage + " Invalid value found in section " + ICSP::section + ".", 
        ICSP::numberOfPointsIntegratedCrossSections + " > 0 must be satisfied."
    );

    // Ensure propagatorIntegralPrecision>0
    bool isPropagatorIntegralPrecisionValid = config.validatePositiveDouble(
        ICSP::section,
        ICSP::propagatorIntegralPrecision,
        invalidFileMessage + " Invalid value found in section " + ICSP::section + ".", 
        ICSP::propagatorIntegralPrecision + " > 0 must be satisfied."
    );

    // Ensure largeAngleScatteringContribution is true or false
    bool isLargeAngleScatteringContribution = config.validateBool(
        ICSP::section,
        ICSP::largeAngleScatteringContribution,
        invalidFileMessage + " Invalid value found in section " + ICSP::section + ".", 
        ICSP::largeAngleScatteringContribution + " must be true (1, yes, on) or false (0, no, off)."
    );

    // Ensure crossSectionIntegralPrecision>0
    bool isCrossSectionIntegralPrecisionValid = config.validatePositiveDouble(
        ICSP::section,
        ICSP::crossSectionIntegralPrecision,
        invalidFileMessage + " Invalid value found in section " + ICSP::section + ".", 
        ICSP::crossSectionIntegralPrecision + " > 0 must be satisfied."
    );

    // Ensure integratedCrossSectionIntegralPrecision_dXdY>0
    bool isIntegratedCrossSectionIntegralPrecision_dXdYValid = config.validatePositiveDouble(
        ICSP::section,
        ICSP::integratedCrossSectionIntegralPrecision_dXdY,
        invalidFileMessage + " Invalid value found in section " + ICSP::section + ".", 
        ICSP::integratedCrossSectionIntegralPrecision_dXdY + " > 0 must be satisfied."
    );

    // Ensure integratedCrossSectionIntegralPrecision_dX>0
    bool isIntegratedCrossSectionIntegralPrecision_dXValid = config.validatePositiveDouble(
        ICSP::section,
        ICSP::integratedCrossSectionIntegralPrecision_dX,
        invalidFileMessage + " Invalid value found in section " + ICSP::section + ".", 
        ICSP::integratedCrossSectionIntegralPrecision_dX + " > 0 must be satisfied."
    );

    // Ensure approximationMethod is valid
    string approximationMethod = config.getValue(ICSP::section, ICSP::approximationMethod);
    bool isApproximationMethodValid = isValidIntegratedCrossSectionApproximationMethod(approximationMethod);
    if( !isApproximationMethodValid ){ cout << invalidFileMessage << endl; }

    // The function return true only if all tests passed
    return isNumberOfPointsIntegratedCrossSectionsValid && 
           isPropagatorIntegralPrecisionValid && 
           isLargeAngleScatteringContribution &&
           isCrossSectionIntegralPrecisionValid &&
           isIntegratedCrossSectionIntegralPrecision_dXdYValid &&
           isIntegratedCrossSectionIntegralPrecision_dXValid &&
           isApproximationMethodValid;
}

bool SU3NJL3DCutoffFileParser::Vacuum::VacuumMasses::validateFile() const
{
    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffConfigKeys::ModelParameters::section, 
        SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section, 
        SU3NJL3DCutoffConfigKeys::VacuumMassesParameters::section
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
    bool areVacuumMassesParametersValid = validateVacuumMassesParameters();

    // The function return true only if all tests passed
    return allRequiredSectionsPresent &&
           areModelParametersValid &&
           areDimensionfulCouplingsValid &&
           areVacuumMassesParametersValid;
}

void SU3NJL3DCutoffFileParser::Vacuum::VacuumMasses::evaluate() const
{
    // Model Parameters
    namespace MP = SU3NJL3DCutoffConfigKeys::ModelParameters;
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);
    
    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMasses
    namespace VMP = SU3NJL3DCutoffConfigKeys::VacuumMassesParameters;
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(
        stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
        cutoff, 
        couplings, 
        upQuarkCurrentMass, 
        downQuarkCurrentMass, 
        strangeQuarkCurrentMass
    );
    parameters.setParameterSetName(parameterSetName);

    SU3NJL3DCutoffVacuum::evaluateVacuumMasses(
        parameters, 
        precisionVacuum, 
        stringToMultiRootFindingMethod(methodVacuum), 
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess
    );
}

bool SU3NJL3DCutoffFileParser::FixedTempRhoBEqualChemPot::FirstOrderLine::validateFile() const
{   
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
    const SU3NJL3DCutoffFileParser::Vacuum::VacuumMasses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffConfigKeys::VacuumToFiniteBaryonDensityParameters::section,
        SU3NJL3DCutoffConfigKeys::FirstOrderLineParameters::section
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
    bool areVacuumToFiniteBaryonDensityParametersValid = validateVacuumToFiniteBaryonDensityParameters();
    bool areFirstOrderLineParametersValid = validateFirstOrderLineParameters();

    // The function return true only if all tests passed
    return vacuumValidations && 
           allRequiredSectionsPresent &&
           areVacuumToFiniteBaryonDensityParametersValid &&
           areFirstOrderLineParametersValid;
}

void SU3NJL3DCutoffFileParser::FixedTempRhoBEqualChemPot::FirstOrderLine::evaluate() const
{
    // Model Parameters
    namespace MP = SU3NJL3DCutoffConfigKeys::ModelParameters;
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);
    
    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMasses
    namespace VMP = SU3NJL3DCutoffConfigKeys::VacuumMassesParameters;
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;

    // VacuumToFiniteBaryonDensity
    namespace VFBDP = SU3NJL3DCutoffConfigKeys::VacuumToFiniteBaryonDensityParameters;
    cout << "\n" << VFBDP::section << ": " << endl;
    
    double minimumBaryonDensity = config.getDouble(VFBDP::section, VFBDP::minimumBaryonDensity);
    double maximumBaryonDensity = config.getDouble(VFBDP::section, VFBDP::maximumBaryonDensity);
    int numberOfPoints = config.getInt(VFBDP::section, VFBDP::numberOfPoints);
    double precisionZeroTempSol = config.getDouble(VFBDP::section, VFBDP::precisionZeroTempSol);
    string methodZeroTempSol = config.getValue(VFBDP::section, VFBDP::methodZeroTempSol);

    cout << VFBDP::minimumBaryonDensity << " = " << minimumBaryonDensity << endl;
    cout << VFBDP::maximumBaryonDensity << " = " << maximumBaryonDensity << endl;
    cout << VFBDP::numberOfPoints << " = " << numberOfPoints << endl;
    cout << VFBDP::precisionZeroTempSol << " = " << precisionZeroTempSol << endl;
    cout << VFBDP::methodZeroTempSol << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodZeroTempSol)) << endl;

    // FirstOrderLine
    namespace FOLP = SU3NJL3DCutoffConfigKeys::FirstOrderLineParameters;
    cout << "\n" << FOLP::section << ":\n";
    
    double precisionTransitionPointSol = config.getDouble(FOLP::section, FOLP::precisionTransitionPointSol);
    string methodTransitionPointSol = config.getValue(FOLP::section, FOLP::methodTransitionPointSol);
    double deltaT = config.getDouble(FOLP::section, FOLP::deltaT);
    double massDifferenceCEP = config.getDouble(FOLP::section, FOLP::massDifferenceCEP);
    
    cout << FOLP::precisionTransitionPointSol  << " = " << precisionTransitionPointSol << endl;
    cout << FOLP::methodTransitionPointSol << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodTransitionPointSol)) << endl;
    cout << FOLP::deltaT << " = " << deltaT << endl;
    cout << FOLP::massDifferenceCEP << " = " << massDifferenceCEP << endl;

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(
        stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
        cutoff, 
        couplings, 
        upQuarkCurrentMass, 
        downQuarkCurrentMass, 
        strangeQuarkCurrentMass
    );
    parameters.setParameterSetName(parameterSetName);

    // Perform calculation
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot::evaluateFirstOrderLine(
        parameters, 
        precisionVacuum, 
        stringToMultiRootFindingMethod(methodVacuum), 
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess, 
        minimumBaryonDensity, 
        maximumBaryonDensity, 
        numberOfPoints, 
        precisionZeroTempSol, 
        stringToMultiRootFindingMethod(methodZeroTempSol), 
        true, 
        precisionTransitionPointSol, 
        stringToMultiRootFindingMethod(methodTransitionPointSol), 
        deltaT, 
        massDifferenceCEP
    );
}

bool SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricCrossSections::validateFile() const
{
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
    const SU3NJL3DCutoffFileParser::Vacuum::VacuumMasses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffConfigKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters::section,
        SU3NJL3DCutoffConfigKeys::FiniteTemperatureToFiniteChemicalPotentialParameters::section,
        SU3NJL3DCutoffConfigKeys::CrossSectionsParameters::section
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
    bool areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid = validateVacuumToFiniteTemperatureAtZeroChemicalPotentialParameters();
    bool areFiniteTemperatureToFiniteChemicalPotentialParametersValid = validateFiniteTemperatureToFiniteChemicalPotentialParameters();
    bool areCrossSectionsParametersValid = validateCrossSectionsParameters();

    // The function return true only if all tests passed
    return vacuumValidations && 
           allRequiredSectionsPresent && 
           areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid && 
           areFiniteTemperatureToFiniteChemicalPotentialParametersValid &&
           areCrossSectionsParametersValid;
}

void SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricCrossSections::evaluate() const
{   
    // Model Parameters
    namespace MP = SU3NJL3DCutoffConfigKeys::ModelParameters;
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);

    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMassesParameters
    namespace VMP = SU3NJL3DCutoffConfigKeys::VacuumMassesParameters;
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;

    // Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(
        stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
        cutoff, 
        couplings, 
        upQuarkCurrentMass, 
        downQuarkCurrentMass, 
        strangeQuarkCurrentMass
    );
    parameters.setParameterSetName(parameterSetName);

    // VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters
    namespace VFTZCPP = SU3NJL3DCutoffConfigKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    cout << "\n" << VFTZCPP::section << ": " << endl;

    double temperature = config.getDouble(VFTZCPP::section, VFTZCPP::temperature);
    int numberOfPointsFromVacToFinTemp = config.getInt(VFTZCPP::section, VFTZCPP::numberOfPointsFromVacToFinTemp);
    double precisionVacToFinTemp = config.getDouble(VFTZCPP::section, VFTZCPP::precisionVacToFinTemp);
    string methodVacToFinTemp = config.getValue(VFTZCPP::section, VFTZCPP::methodVacToFinTemp);

    cout << VFTZCPP::temperature << " = " << temperature << endl;
    cout << VFTZCPP::numberOfPointsFromVacToFinTemp << " = " << numberOfPointsFromVacToFinTemp << endl;
    cout << VFTZCPP::precisionVacToFinTemp << " = " << precisionVacToFinTemp << endl;
    cout << VFTZCPP::methodVacToFinTemp << " = " << methodVacToFinTemp << endl;

    // FiniteTemperatureToFiniteChemicalPotentialParameters
    namespace FTFCPP = SU3NJL3DCutoffConfigKeys::FiniteTemperatureToFiniteChemicalPotentialParameters;
    cout << "\n" << FTFCPP::section << ": " << endl;
    
    double chemPot = config.getDouble(FTFCPP::section, FTFCPP::chemPot);
    int numberOfPointsFromFinTempToFinChemPot = config.getInt(FTFCPP::section, FTFCPP::numberOfPointsFromFinTempToFinChemPot);
    double precisionFinTempToFinChemPot = config.getDouble(FTFCPP::section, FTFCPP::precisionFinTempToFinChemPot);
    string methodFinTempToFinChemPot = config.getValue(FTFCPP::section, FTFCPP::methodFinTempToFinChemPot);

    cout << FTFCPP::chemPot << " = " << chemPot << endl;
    cout << FTFCPP::numberOfPointsFromFinTempToFinChemPot << " = " << numberOfPointsFromFinTempToFinChemPot << endl;
    cout << FTFCPP::precisionFinTempToFinChemPot << " = " << precisionFinTempToFinChemPot << endl;
    cout << FTFCPP::methodFinTempToFinChemPot << " = " << methodFinTempToFinChemPot << endl;

    // CrossSectionsParameters
    namespace CSP = SU3NJL3DCutoffConfigKeys::CrossSectionsParameters;
    cout << "\n" << CSP::section << ": " << endl;

    double propagatorIntegralPrecision = config.getDouble(CSP::section, CSP::propagatorIntegralPrecision);
    bool largeAngleScatteringContribution = config.getBool(CSP::section, CSP::largeAngleScatteringContribution);
    double precisionCrossSections = config.getDouble(CSP::section, CSP::precisionCrossSections);
    int numberOfPointsCrossSections = config.getInt(CSP::section, CSP::numberOfPointsCrossSections);
    int numberOfThreads = config.getInt(CSP::section, CSP::numberOfThreads);

    cout << CSP::propagatorIntegralPrecision << " = " << propagatorIntegralPrecision << endl;
    cout << CSP::largeAngleScatteringContribution << " = " << largeAngleScatteringContribution << endl;
    cout << CSP::precisionCrossSections << " = " << precisionCrossSections << endl;
    cout << CSP::numberOfPointsCrossSections << " = " << numberOfPointsCrossSections << endl;
    cout << CSP::numberOfThreads << " = " << numberOfThreads << endl;

    // Calculate Cross Sections
    SU3NJL3DCutoffFixedChemPotTemp::evaluateCrossSectionsEqualLightMasses(
    	parameters, 
		precisionVacuum, 
		stringToMultiRootFindingMethod(methodVacuum), 
		upQuarkMassGuess, 
		strangeQuarkMassGuess,
        temperature,
        numberOfPointsFromVacToFinTemp,
        precisionVacToFinTemp,
        stringToMultiRootFindingMethod(methodVacToFinTemp),
        chemPot,
        numberOfPointsFromFinTempToFinChemPot,
        precisionFinTempToFinChemPot,
        stringToMultiRootFindingMethod(methodFinTempToFinChemPot),
        propagatorIntegralPrecision,
        largeAngleScatteringContribution,
        precisionCrossSections,
        numberOfPointsCrossSections,
        numberOfThreads
    );
}

bool SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricIntegratedCrossSections::validateFileZeroChemicalPotential() const
{
    namespace VFT = SU3NJL3DCutoffConfigKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    namespace LHT = SU3NJL3DCutoffConfigKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters;
    namespace ICS = SU3NJL3DCutoffConfigKeys::IntegratedCrossSectionsParameters;

    bool areTemperaturesEqual = (config.getValue(VFT::section, VFT::temperature)==config.getValue(LHT::section, LHT::minimumTemp));
    if ( !areTemperaturesEqual )
    {
        cout << invalidFileMessage + " Invalid values found in sections " + VFT::section + " and " + LHT::section + ": \n";
        cout << VFT::temperature << "==" << LHT::minimumTemp << " must be satisfied.\n";
    }

    bool areNumberOfPointsEqual = (config.getValue(LHT::section, LHT::numberOfPointsFromLowToHighTemp)==config.getValue(ICS::section, ICS::numberOfPointsIntegratedCrossSections));
    if ( !areNumberOfPointsEqual )
    {
        cout << invalidFileMessage + " Invalid values found in sections " + LHT::section + " and " + ICS::section + ": \n";
        cout << LHT::numberOfPointsFromLowToHighTemp << "==" << ICS::numberOfPointsIntegratedCrossSections << " must be satisfied.\n";
    }
    
    return areTemperaturesEqual && areNumberOfPointsEqual;
}

bool SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricIntegratedCrossSections::validateFile(string physicalScenario) const
{   
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
	const SU3NJL3DCutoffFileParser::Vacuum::VacuumMasses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffConfigKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters::section,
        SU3NJL3DCutoffConfigKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters::section,
        SU3NJL3DCutoffConfigKeys::IntegratedCrossSectionsParameters::section
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
    bool areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid = validateVacuumToFiniteTemperatureAtZeroChemicalPotentialParameters();
    bool areLowToHighTemperatureAtZeroChemicalPotentialParametersValid = validateLowToHighTemperatureAtZeroChemicalPotentialParameters();
    bool areIntegratedCrossSectionsParameters = validateIntegratedCrossSectionsParameters();

    // Physical scenario specific validations
    bool isPhysicalScenarioValid = true;
    if ( physicalScenario==zeroChemicalPotential )
    {
        isPhysicalScenarioValid = validateFileZeroChemicalPotential();
    }

    return vacuumValidations && 
           allRequiredSectionsPresent &&
           areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid &&
           areLowToHighTemperatureAtZeroChemicalPotentialParametersValid &&
           areIntegratedCrossSectionsParameters &&
           isPhysicalScenarioValid;
}

void SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricIntegratedCrossSections::evaluate(string physicalScenario) const
{   
    namespace MP = SU3NJL3DCutoffConfigKeys::ModelParameters;
    namespace VMP = SU3NJL3DCutoffConfigKeys::VacuumMassesParameters;
    namespace VFT = SU3NJL3DCutoffConfigKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    namespace LHT = SU3NJL3DCutoffConfigKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters;
    namespace ICS = SU3NJL3DCutoffConfigKeys::IntegratedCrossSectionsParameters;

    // Model Parameters
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);

    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffConfigKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMassesParameters
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;
    
    // VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters
    cout << "\n" << VFT::section << ": " << endl;

    double temperature = config.getDouble(VFT::section, VFT::temperature);
    int numberOfPointsFromVacToFinTemp = config.getInt(VFT::section, VFT::numberOfPointsFromVacToFinTemp);
    double precisionVacToFinTemp = config.getDouble(VFT::section, VFT::precisionVacToFinTemp);
    string methodVacToFinTemp = config.getValue(VFT::section, VFT::methodVacToFinTemp);

    cout << VFT::temperature << " = " << temperature << endl;
    cout << VFT::numberOfPointsFromVacToFinTemp << " = " << numberOfPointsFromVacToFinTemp << endl;
    cout << VFT::precisionVacToFinTemp << " = " << precisionVacToFinTemp << endl;
    cout << VFT::methodVacToFinTemp << " = " << methodVacToFinTemp << endl;

    // LowToHighTemperatureAtZeroChemicalPotentialParameters
    cout << "\n" << LHT::section << ": " << endl;
    double minimumTemp = config.getDouble(LHT::section, LHT::minimumTemp);
    double maximumTemp = config.getDouble(LHT::section, LHT::maximumTemp);
    double numberOfPointsFromLowToHighTemp = config.getInt(LHT::section, LHT::numberOfPointsFromLowToHighTemp);
    double precisionLowToHighTemp = config.getDouble(LHT::section, LHT::precisionLowToHighTemp);
    string methodLowToHighTemp = config.getValue(LHT::section, LHT::methodLowToHighTemp);  

    cout << LHT::minimumTemp << " = " << minimumTemp << endl;
    cout << LHT::maximumTemp << " = " << maximumTemp << endl;
    cout << LHT::numberOfPointsFromLowToHighTemp << " = " << numberOfPointsFromLowToHighTemp << endl;
    cout << LHT::precisionLowToHighTemp << " = " << precisionLowToHighTemp << endl;
    cout << LHT::methodLowToHighTemp << " = " << methodLowToHighTemp << endl;

    // IntegratedCrossSectionsParameters
    cout << "\n" << ICS::section << ": " << endl;
    int numberOfPointsIntegratedCrossSections = config.getInt(ICS::section, ICS::numberOfPointsIntegratedCrossSections);
    double propagatorIntegralPrecision = config.getDouble(ICS::section, ICS::propagatorIntegralPrecision);
    bool largeAngleScatteringContribution = config.getBool(ICS::section, ICS::largeAngleScatteringContribution);
    double crossSectionIntegralPrecision = config.getDouble(ICS::section, ICS::crossSectionIntegralPrecision);
    double integratedCrossSectionIntegralPrecision_dXdY = config.getDouble(ICS::section, ICS::integratedCrossSectionIntegralPrecision_dXdY);
    double integratedCrossSectionIntegralPrecision_dX = config.getDouble(ICS::section, ICS::integratedCrossSectionIntegralPrecision_dX);
    string approximationMethod = config.getValue(ICS::section, ICS::approximationMethod);
    //int numberOfThreads = config.getInt(ICS::section, SU3NJL3DCutoffConfigKeys::CrossSectionsParameters::numberOfThreads);
  
    cout << ICS::numberOfPointsIntegratedCrossSections << " = " << numberOfPointsIntegratedCrossSections << endl;
    cout << ICS::propagatorIntegralPrecision << " = " << propagatorIntegralPrecision << endl;
    cout << ICS::largeAngleScatteringContribution << " = " << largeAngleScatteringContribution << endl;
    cout << ICS::crossSectionIntegralPrecision << " = " << crossSectionIntegralPrecision << endl;
    cout << ICS::integratedCrossSectionIntegralPrecision_dXdY << " = " << integratedCrossSectionIntegralPrecision_dXdY << endl;
    cout << ICS::integratedCrossSectionIntegralPrecision_dX << " = " << integratedCrossSectionIntegralPrecision_dX << endl;
    cout << ICS::approximationMethod << " = " << approximationMethod << endl;
    //cout << SU3NJL3DCutoffConfigKeys::CrossSectionsParameters::numberOfThreads << " = " << numberOfThreads << endl;

    // Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(
        stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
        cutoff, 
        couplings, 
        upQuarkCurrentMass, 
        downQuarkCurrentMass, 
        strangeQuarkCurrentMass
    );
    parameters.setParameterSetName(parameterSetName);
    
    if ( physicalScenario==zeroChemicalPotential )
    {
        evaluateIsospinSymmetricIntegratedCrossSectionsWithZeroChemicalPotential(
            parameters, 
            precisionVacuum, 
            stringToMultiRootFindingMethod(methodVacuum), 
            upQuarkMassGuess, 
            strangeQuarkMassGuess, 
            precisionVacToFinTemp,
            stringToMultiRootFindingMethod(methodVacToFinTemp), 
            minimumTemp, 
            maximumTemp, 
            numberOfPointsFromVacToFinTemp, 
            numberOfPointsFromLowToHighTemp, 
            largeAngleScatteringContribution, 
            stringToIntegratedCrossSectionApproximationMethod(approximationMethod),
            propagatorIntegralPrecision,
            crossSectionIntegralPrecision,
            integratedCrossSectionIntegralPrecision_dXdY,
            integratedCrossSectionIntegralPrecision_dX
        );
    }
}

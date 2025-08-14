#include <iostream>
#include <cmath>
#include "gsl_wrapper/root_solver_gsl.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "njl_model/njl_regularization_schemes.h"
#include "njl_model/NJLDimensionlessCouplings.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedTempRhoBEqualChemPot.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"

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

    // Ensure chemPot>0
    bool isChemPotValid = config.validatePositiveDouble(
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

    // The function return true only if all tests passed
    return isPropagatorIntegralPrecisionValid && 
           isPrecisionCrossSectionsValid &&
           isNumberOfPointsCrossSectionsValid;
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

    cout << CSPKeys::propagatorIntegralPrecision << " = " << propagatorIntegralPrecision << endl;
    cout << CSPKeys::largeAngleScatteringContribution << " = " << largeAngleScatteringContribution << endl;
    cout << CSPKeys::precisionCrossSections << " = " << precisionCrossSections << endl;
    cout << CSPKeys::numberOfPointsCrossSections << " = " << numberOfPointsCrossSections << endl;

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
        numberOfPointsCrossSections
    );
}

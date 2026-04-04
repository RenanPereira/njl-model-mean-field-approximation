#include <iostream>
#include <cmath>

#include "gsl_wrapper/root_solver_gsl.h"

#include "njl_model/njl_regularization_schemes.h"
#include "njl_model/NJLDimensionlessCouplings.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParserKeys.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedTempRhoBEqualChemPot.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.h"

using namespace std;

namespace SU3NJL3DCutoffFileParser
{

NJLDimensionfulCouplings Common::extractDimensionfulCouplings() const
{	
	double cutoff = config.getDouble(
        SU3NJL3DCutoffFileParserKeys::ModelParameters::section,                                   
        SU3NJL3DCutoffFileParserKeys::ModelParameters::cutoff
    );

    LagrangianInteractions interaction = stringToLagrangianInteractions(
        config.getValue(
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section,             
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::lagrangianInteractions
        )
    );
	cout << "interaction = " << toString(interaction) << endl;

	if ( interaction==SP4Q_DET2NFQ )
	{
		double fourQuarkSPCouplingCutoff2 = config.getDouble(
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section,              
            NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING
        );
		double determinantCouplingCutoff5 = config.getDouble(
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section,              
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
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::FOUR_QUARK_SP_COUPLING
        );
		double determinantCouplingCutoff5 = config.getDouble(
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::DETERMINANT_COUPLING
        );
		double eightQuarkSPOziViolatingCouplingCutoff8 = config.getDouble(
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section, 
            NJLDimensionlessCouplings::EIGHT_QUARK_SP_OZI_COUPLING
        );
		double eightQuarkSPNonOziViolatingCouplingCutoff8 = config.getDouble(
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section, 
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

bool Common::validateModelParameters() const
{   
    // Validate regularizationScheme
    string regularizationScheme = config.getValue(
        SU3NJL3DCutoffFileParserKeys::ModelParameters::section, 
        SU3NJL3DCutoffFileParserKeys::ModelParameters::regularizationScheme
    );
    bool isRegularizationSchemeValid = isValidNJL3DCutoffRegularizationScheme(regularizationScheme);
    if( !isRegularizationSchemeValid ){ cout << invalidFileMessage << endl; }

    // Ensure cutoff>0
    bool isCutoffValid = config.validatePositiveDouble(
        SU3NJL3DCutoffFileParserKeys::ModelParameters::section, 
        SU3NJL3DCutoffFileParserKeys::ModelParameters::cutoff, 
        invalidFileMessage, 
        SU3NJL3DCutoffFileParserKeys::ModelParameters::cutoff + " > 0 must be satisfied."
    );

    // Ensure upQuarkCurrentMass>0
    bool isUpQuarkCurrentMassValid = config.validatePositiveDouble(
        SU3NJL3DCutoffFileParserKeys::ModelParameters::section, 
        SU3NJL3DCutoffFileParserKeys::ModelParameters::upQuarkCurrentMass, 
        invalidFileMessage, 
        SU3NJL3DCutoffFileParserKeys::ModelParameters::upQuarkCurrentMass + " > 0 must be satisfied."
    );

    // Ensure downQuarkCurrentMassInGeV>0
    bool isDownQuarkCurrentMassValid = config.validatePositiveDouble(
        SU3NJL3DCutoffFileParserKeys::ModelParameters::section, 
        SU3NJL3DCutoffFileParserKeys::ModelParameters::downQuarkCurrentMass, 
        invalidFileMessage, 
        SU3NJL3DCutoffFileParserKeys::ModelParameters::downQuarkCurrentMass + " > 0 must be satisfied."
    );

    // Ensure strangeQuarkCurrentMass>0
    bool isStrangeQuarkCurrentMassValid = config.validatePositiveDouble(
        SU3NJL3DCutoffFileParserKeys::ModelParameters::section, 
        SU3NJL3DCutoffFileParserKeys::ModelParameters::strangeQuarkCurrentMass, 
        invalidFileMessage, 
        SU3NJL3DCutoffFileParserKeys::ModelParameters::strangeQuarkCurrentMass + " > 0 must be satisfied."
    );

    // The function return true only if all tests passed
    return isRegularizationSchemeValid &&
           isCutoffValid &&
           isUpQuarkCurrentMassValid &&
           isDownQuarkCurrentMassValid &&
           isStrangeQuarkCurrentMassValid;
}

bool Common::validateDimensionfulCouplings() const
{
    // Validate lagrangianInteractions
    string lagrangianInteractions = config.getValue(
        SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section,                                             
        SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::lagrangianInteractions
    );
    bool isLagrangianInteractionsValid = isValidLagrangianInteractions(lagrangianInteractions);
    if( !isLagrangianInteractionsValid ){ cout << invalidFileMessage << endl; }

    bool areCouplingsValid = false;
    if ( isLagrangianInteractionsValid )
    {
        areCouplingsValid = validateNJLDimensionfulCouplings(
            config,
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section,
            SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::lagrangianInteractions
        );
    }

    return areCouplingsValid;
}

bool Common::validateVacuumMassesParameters() const
{   
    // Ensure precisionVacuum>0
    bool isPrecisionVacuumValid = config.validatePositiveDouble(
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::section,
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::precisionVacuum, 
        invalidFileMessage, 
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::precisionVacuum + " > 0 must be satisfied."
    );

    // Ensure methodVacuum is valid
    string methodVacuum = config.getValue(
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::methodVacuum
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodVacuum);
    if( !isRootFindingMethodValid ){ cout << invalidFileMessage << endl; }

    // Ensure up quark mass is valid
    bool isUpQuarkMassGuessValid = config.validatePositiveDouble(
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::upQuarkMassGuess, 
        invalidFileMessage, 
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::upQuarkMassGuess + " > 0 must be satisfied."
    );

    // Ensure down quark mass is valid
    bool isDownQuarkMassGuessValid = config.validatePositiveDouble(
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::downQuarkMassGuess, 
        invalidFileMessage, 
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::downQuarkMassGuess + " > 0 must be satisfied."
    );

    // Ensure strange quark mass is valid
    bool isStrangeQuarkMassGuessValid = config.validatePositiveDouble(
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::section, 
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::strangeQuarkMassGuess, 
        invalidFileMessage, 
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::strangeQuarkMassGuess + " > 0 must be satisfied."
    );

    // The function returns true only if all tests passed
    return isPrecisionVacuumValid &&
           isRootFindingMethodValid &&
           isDownQuarkMassGuessValid &&
           isStrangeQuarkMassGuessValid &&
           isUpQuarkMassGuessValid;
}

bool Common::validateVacuumToFiniteBaryonDensityParameters() const
{   
    namespace VFBD = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteBaryonDensityParameters;

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

bool Common::validateFirstOrderLineParameters() const
{   
    namespace FOLP = SU3NJL3DCutoffFileParserKeys::FirstOrderLineParameters;

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

bool Common::validateVacuumToFiniteChemicalPotentialParameters() const
{
    namespace VTFCPP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteChemicalPotentialParameters;
    const string sectionError = invalidFileMessage + " Invalid value found in section " + VTFCPP::section + ".";

    // Ensure chemicalPotential>0
    bool isChemicalPotentialValid = config.validateNonNegativeDouble(
        VTFCPP::section,
        VTFCPP::chemicalPotential,
        sectionError, 
        VTFCPP::chemicalPotential + " >= 0 must be satisfied."
    );

    // Ensure numberOfPointsVacToChemPot>0
    bool isNumberOfPointsVacToChemPotValid = config.validatePositiveInteger(
        VTFCPP::section,
        VTFCPP::numberOfPointsVacToChemPot, 
        sectionError, 
        VTFCPP::numberOfPointsVacToChemPot + " > 0 must be satisfied."
    );

    // Ensure precisionVacToFinTemp>0
    bool isPrecisionVacToChemPotValid = config.validatePositiveDouble(
        VTFCPP::section,
        VTFCPP::precisionVacToChemPot,
        sectionError, 
        VTFCPP::precisionVacToChemPot + " > 0 must be satisfied."
    );

    // Ensure methodVacToChemPot is valid
    string methodVacToChemPot = config.getValue(
        VTFCPP::section, 
        VTFCPP::methodVacToChemPot
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodVacToChemPot, sectionError);

    // The function return true only if all tests passed
    return isChemicalPotentialValid && 
           isNumberOfPointsVacToChemPotValid &&
           isPrecisionVacToChemPotValid &&
           isRootFindingMethodValid;
}

bool Common::validateVacuumToTemperatureParameters() const
{
    namespace VTTP = SU3NJL3DCutoffFileParserKeys::VacuumToTemperatureParameters;
    const string sectionError = invalidFileMessage + " Invalid value found in section " + VTTP::section + ".";

    // Ensure temperature>0
    bool isTemperatureValid = config.validateNonNegativeDouble(
        VTTP::section,
        VTTP::temperature,
        sectionError, 
        VTTP::temperature + " >= 0 must be satisfied."
    );

    // Ensure numberOfPointsVacToTemp>0
    bool isNumberOfPointsVacToTemp = config.validatePositiveInteger(
        VTTP::section,
        VTTP::numberOfPointsVacToTemp, 
        sectionError, 
        VTTP::numberOfPointsVacToTemp + " > 0 must be satisfied."
    );

    // Ensure precisionVacToTemp>0
    bool isPrecisionVacToTempValid = config.validatePositiveDouble(
        VTTP::section,
        VTTP::precisionVacToTemp,
        sectionError, 
        VTTP::precisionVacToTemp + " > 0 must be satisfied."
    );

    // Ensure methodVacToTemp is valid
    string methodVacToTemp = config.getValue(
        VTTP::section, 
        VTTP::methodVacToTemp
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodVacToTemp, sectionError);

    // The function return true only if all tests passed
    return isTemperatureValid && 
           isNumberOfPointsVacToTemp &&
           isPrecisionVacToTempValid &&
           isRootFindingMethodValid;
}

bool Common::validateToTemperatureParameters() const
{
    namespace TTP = SU3NJL3DCutoffFileParserKeys::ToTemperatureParameters;
    const string sectionError = invalidFileMessage + " Invalid value found in section " + TTP::section + ".";

    // Ensure temperature>0
    bool isTemperatureValid = config.validatePositiveDouble(
        TTP::section,
        TTP::temperature,
        sectionError, 
        TTP::temperature + " > 0 must be satisfied."
    );

    // Ensure numberOfPointsUpToTemp>0
    bool isNumberOfPointsUpToTempValid = config.validatePositiveInteger(
        TTP::section,
        TTP::numberOfPointsUpToTemp, 
        sectionError, 
        TTP::numberOfPointsUpToTemp + " > 0 must be satisfied."
    );

    // Ensure precisionUpToTemp>0
    bool isPrecisionUpToTempValid = config.validatePositiveDouble(
        TTP::section,
        TTP::precisionUpToTemp,
        sectionError, 
        TTP::precisionUpToTemp + " > 0 must be satisfied."
    );

    // Ensure methodUpToTemp is valid
    string methodUpToTemp = config.getValue(
        TTP::section, 
        TTP::methodUpToTemp
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodUpToTemp, sectionError);

    // The function return true only if all tests passed
    return isTemperatureValid && 
           isNumberOfPointsUpToTempValid &&
           isPrecisionUpToTempValid &&
           isRootFindingMethodValid;
}

bool Common::validateToChemicalPotentialSymmetricParameters() const
{
    namespace TCPSP = SU3NJL3DCutoffFileParserKeys::ToChemicalPotentialSymmetricParameters;
    const string sectionError = invalidFileMessage + " Invalid value found in section " + TCPSP::section + ".";

    // Ensure chemicalPotential>0
    bool isChemicalPotentialValid = config.validatePositiveDouble(
        TCPSP::section,
        TCPSP::chemicalPotential,
        sectionError, 
        TCPSP::chemicalPotential + " > 0 must be satisfied."
    );

    // Ensure numberOfPointsToChemPot>0
    bool isNumberOfPointsToChemPotValid = config.validatePositiveInteger(
        TCPSP::section,
        TCPSP::numberOfPointsToChemPot, 
        sectionError, 
        TCPSP::numberOfPointsToChemPot + " > 0 must be satisfied."
    );

    // Ensure precisionToChemPot>0
    bool isPrecisionToChemPotValid = config.validatePositiveDouble(
        TCPSP::section,
        TCPSP::precisionToChemPot,
        sectionError, 
        TCPSP::precisionToChemPot + " > 0 must be satisfied."
    );

    // Ensure methodToChemPot is valid
    string methodToChemPot = config.getValue(
        TCPSP::section, 
        TCPSP::methodToChemPot
    );
    bool isRootFindingMethodValid = isValidMultiRootFindingMethod(methodToChemPot, sectionError);

    // The function return true only if all tests passed
    return isChemicalPotentialValid && 
           isNumberOfPointsToChemPotValid &&
           isPrecisionToChemPotValid &&
           isRootFindingMethodValid;
}

bool Common::validateVacuumToFiniteTemperatureAtZeroChemicalPotentialParameters() const
{
    namespace VFTZCPP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;

    // Ensure nearVacuumTemperature>=0
    bool isNearVacuumTemperatureValid = config.validateNonNegativeDouble(
        VFTZCPP::section,
        VFTZCPP::nearVacuumTemperature ,
        invalidFileMessage + " Invalid value found in section " + VFTZCPP::section + ".", 
        VFTZCPP::nearVacuumTemperature + " >= 0 must be satisfied."
    );

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
    return isNearVacuumTemperatureValid &&
           isTemperatureValid && 
           isNumberOfPointsFromVacToFinTempValid &&
           isPrecisionVacToFinTempValid &&
           isRootFindingMethodValid;
}

bool Common::validateLowToHighTemperatureAtZeroChemicalPotentialParameters() const
{   
    namespace LHTZCPP = SU3NJL3DCutoffFileParserKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters;

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

bool Common::validateFiniteTemperatureToFiniteChemicalPotentialParameters() const
{
    namespace FTFCPP = SU3NJL3DCutoffFileParserKeys::FiniteTemperatureToFiniteChemicalPotentialParameters;

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

bool Common::validateCrossSectionsParameters() const
{
    namespace CSP = SU3NJL3DCutoffFileParserKeys::CrossSectionsParameters;

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

bool Common::validateIntegratedCrossSectionsParameters() const
{
    namespace ICSP = SU3NJL3DCutoffFileParserKeys::IntegratedCrossSectionsParameters;

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

    // Ensure numberOfThreads>0
    bool isNumberOfThreadsValid = config.validatePositiveInteger(
        ICSP::section,
        ICSP::numberOfThreads, 
        invalidFileMessage + " Invalid value found in section " + ICSP::section + ".", 
        ICSP::numberOfThreads + " > 0 must be satisfied."
    );

    // The function return true only if all tests passed
    return isNumberOfPointsIntegratedCrossSectionsValid && 
           isPropagatorIntegralPrecisionValid && 
           isLargeAngleScatteringContribution &&
           isCrossSectionIntegralPrecisionValid &&
           isIntegratedCrossSectionIntegralPrecision_dXdYValid &&
           isIntegratedCrossSectionIntegralPrecision_dXValid &&
           isApproximationMethodValid &&
           isNumberOfThreadsValid;
}

bool Common::checkRequiredSections(const vector<string> requiredSections) const
{
    // Check for missing sections
    bool allRequiredSectionsPresent = true;

    for (int i = 0; i < int(requiredSections.size()); ++i) 
    {
        string section = requiredSections[i];
        if (config.getSectionsData(section).empty()) 
        {   
            allRequiredSectionsPresent = false;
            cout << "Missing required section: " << section << endl;
        }
    }

    return allRequiredSectionsPresent;
}

void Common::printQualityCheckFailedMessage() const
{
    cout << "The quality check failed for the file: " << config.getFilename() << endl;
}

}

namespace SU3NJL3DCutoffFileParser::Vacuum
{

bool Masses::validateFile() const
{
    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffFileParserKeys::ModelParameters::section, 
        SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section, 
        SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters::section
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

void Masses::evaluate() const
{
    // Model Parameters
    namespace MP = SU3NJL3DCutoffFileParserKeys::ModelParameters;
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);
    
    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toString(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMasses
    namespace VMP = SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters;
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toString(stringToMultiRootFindingMethod(methodVacuum)) << endl;
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

}

namespace SU3NJL3DCutoffFileParser::FixedTempRhoBEqualChemPot
{

bool FirstOrderLine::validateFile() const
{   
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
    const SU3NJL3DCutoffFileParser::Vacuum::Masses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffFileParserKeys::VacuumToFiniteBaryonDensityParameters::section,
        SU3NJL3DCutoffFileParserKeys::FirstOrderLineParameters::section
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

void FirstOrderLine::evaluate() const
{
    // Model Parameters
    namespace MP = SU3NJL3DCutoffFileParserKeys::ModelParameters;
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);
    
    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toString(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMasses
    namespace VMP = SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters;
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toString(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;

    // VacuumToFiniteBaryonDensity
    namespace VFBDP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteBaryonDensityParameters;
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
    cout << VFBDP::methodZeroTempSol << " = " << toString(stringToMultiRootFindingMethod(methodZeroTempSol)) << endl;

    // FirstOrderLine
    namespace FOLP = SU3NJL3DCutoffFileParserKeys::FirstOrderLineParameters;
    cout << "\n" << FOLP::section << ":\n";
    
    double precisionTransitionPointSol = config.getDouble(FOLP::section, FOLP::precisionTransitionPointSol);
    string methodTransitionPointSol = config.getValue(FOLP::section, FOLP::methodTransitionPointSol);
    double deltaT = config.getDouble(FOLP::section, FOLP::deltaT);
    double massDifferenceCEP = config.getDouble(FOLP::section, FOLP::massDifferenceCEP);
    
    cout << FOLP::precisionTransitionPointSol  << " = " << precisionTransitionPointSol << endl;
    cout << FOLP::methodTransitionPointSol << " = " << toString(stringToMultiRootFindingMethod(methodTransitionPointSol)) << endl;
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

}
namespace SU3NJL3DCutoffFileParser::FixedChemPotTemp 
{

bool IsospinSymmetricCrossSections::validateFile() const
{
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
    const SU3NJL3DCutoffFileParser::Vacuum::Masses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters::section,
        SU3NJL3DCutoffFileParserKeys::FiniteTemperatureToFiniteChemicalPotentialParameters::section,
        SU3NJL3DCutoffFileParserKeys::CrossSectionsParameters::section
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

void IsospinSymmetricCrossSections::evaluate() const
{   
    // Model Parameters
    namespace MP = SU3NJL3DCutoffFileParserKeys::ModelParameters;
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);

    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toString(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMassesParameters
    namespace VMP = SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters;
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toString(stringToMultiRootFindingMethod(methodVacuum)) << endl;
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
    namespace VFTZCPP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    cout << "\n" << VFTZCPP::section << ": " << endl;

    double nearVacuumTemperature = config.getDouble(VFTZCPP::section, VFTZCPP::nearVacuumTemperature);
    double temperature = config.getDouble(VFTZCPP::section, VFTZCPP::temperature);
    int numberOfPointsFromVacToFinTemp = config.getInt(VFTZCPP::section, VFTZCPP::numberOfPointsFromVacToFinTemp);
    double precisionVacToFinTemp = config.getDouble(VFTZCPP::section, VFTZCPP::precisionVacToFinTemp);
    string methodVacToFinTemp = config.getValue(VFTZCPP::section, VFTZCPP::methodVacToFinTemp);

    cout << VFTZCPP::nearVacuumTemperature << " = " << nearVacuumTemperature << endl;
    cout << VFTZCPP::temperature << " = " << temperature << endl;
    cout << VFTZCPP::numberOfPointsFromVacToFinTemp << " = " << numberOfPointsFromVacToFinTemp << endl;
    cout << VFTZCPP::precisionVacToFinTemp << " = " << precisionVacToFinTemp << endl;
    cout << VFTZCPP::methodVacToFinTemp << " = " << methodVacToFinTemp << endl;

    // FiniteTemperatureToFiniteChemicalPotentialParameters
    namespace FTFCPP = SU3NJL3DCutoffFileParserKeys::FiniteTemperatureToFiniteChemicalPotentialParameters;
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
    namespace CSP = SU3NJL3DCutoffFileParserKeys::CrossSectionsParameters;
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
    SU3NJL3DCutoffFixedChemPotTemp::evaluateIsospinSymmetricCrossSections(
    	parameters, 
		precisionVacuum, 
		stringToMultiRootFindingMethod(methodVacuum), 
		upQuarkMassGuess, 
		strangeQuarkMassGuess,
        nearVacuumTemperature,
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

bool IsospinSymmetricIntegratedCrossSectionsZeroChemPot::validateTemperatureAndGridConsistency() const
{
    namespace VFTZCPP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    namespace LHT = SU3NJL3DCutoffFileParserKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters;
    namespace ICS = SU3NJL3DCutoffFileParserKeys::IntegratedCrossSectionsParameters;

    bool areTemperaturesEqual = (config.getValue(VFTZCPP::section, VFTZCPP::temperature)==config.getValue(LHT::section, LHT::minimumTemp));
    if ( !areTemperaturesEqual )
    {
        cout << invalidFileMessage + " Invalid values found in sections " + VFTZCPP::section + " and " + LHT::section + ": \n";
        cout << VFTZCPP::temperature << "==" << LHT::minimumTemp << " must be satisfied.\n";
    }

    bool areNumberOfPointsEqual = (config.getValue(LHT::section, LHT::numberOfPointsFromLowToHighTemp)==config.getValue(ICS::section, ICS::numberOfPointsIntegratedCrossSections));
    if ( !areNumberOfPointsEqual )
    {
        cout << invalidFileMessage + " Invalid values found in sections " + LHT::section + " and " + ICS::section + ": \n";
        cout << LHT::numberOfPointsFromLowToHighTemp << "==" << ICS::numberOfPointsIntegratedCrossSections << " must be satisfied.\n";
    }
    
    return areTemperaturesEqual && areNumberOfPointsEqual;
}

bool IsospinSymmetricIntegratedCrossSectionsZeroChemPot::validateFile() const
{   
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
	const SU3NJL3DCutoffFileParser::Vacuum::Masses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters::section,
        SU3NJL3DCutoffFileParserKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters::section,
        SU3NJL3DCutoffFileParserKeys::IntegratedCrossSectionsParameters::section
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
    bool areIntegratedCrossSectionsParametersValid = validateIntegratedCrossSectionsParameters();

    // Other validations
    bool areTemperatureAndGridValid = validateTemperatureAndGridConsistency();

    return vacuumValidations && 
           allRequiredSectionsPresent &&
           areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid &&
           areLowToHighTemperatureAtZeroChemicalPotentialParametersValid &&
           areIntegratedCrossSectionsParametersValid &&
           areTemperatureAndGridValid;
}

void IsospinSymmetricIntegratedCrossSectionsZeroChemPot::evaluate() const
{   
    namespace MP = SU3NJL3DCutoffFileParserKeys::ModelParameters;
    namespace VMP = SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters;
    namespace VFTZCPP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    namespace LHT = SU3NJL3DCutoffFileParserKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters;
    namespace ICSP = SU3NJL3DCutoffFileParserKeys::IntegratedCrossSectionsParameters;

    // Model Parameters
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);

    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toString(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMassesParameters
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toString(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;
    
    // VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters
    cout << "\n" << VFTZCPP::section << ": " << endl;

    double nearVacuumTemperature = config.getDouble(VFTZCPP::section, VFTZCPP::nearVacuumTemperature);
    double temperature = config.getDouble(VFTZCPP::section, VFTZCPP::temperature);
    int numberOfPointsFromVacToFinTemp = config.getInt(VFTZCPP::section, VFTZCPP::numberOfPointsFromVacToFinTemp);
    double precisionVacToFinTemp = config.getDouble(VFTZCPP::section, VFTZCPP::precisionVacToFinTemp);
    string methodVacToFinTemp = config.getValue(VFTZCPP::section, VFTZCPP::methodVacToFinTemp);

    cout << VFTZCPP::nearVacuumTemperature << " = " << nearVacuumTemperature << endl;
    cout << VFTZCPP::temperature << " = " << temperature << endl;
    cout << VFTZCPP::numberOfPointsFromVacToFinTemp << " = " << numberOfPointsFromVacToFinTemp << endl;
    cout << VFTZCPP::precisionVacToFinTemp << " = " << precisionVacToFinTemp << endl;
    cout << VFTZCPP::methodVacToFinTemp << " = " << methodVacToFinTemp << endl;

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
    cout << "\n" << ICSP::section << ": " << endl;
    int numberOfPointsIntegratedCrossSections = config.getInt(ICSP::section, ICSP::numberOfPointsIntegratedCrossSections);
    double propagatorIntegralPrecision = config.getDouble(ICSP::section, ICSP::propagatorIntegralPrecision);
    bool largeAngleScatteringContribution = config.getBool(ICSP::section, ICSP::largeAngleScatteringContribution);
    double crossSectionIntegralPrecision = config.getDouble(ICSP::section, ICSP::crossSectionIntegralPrecision);
    double integratedCrossSectionIntegralPrecision_dXdY = config.getDouble(ICSP::section, ICSP::integratedCrossSectionIntegralPrecision_dXdY);
    double integratedCrossSectionIntegralPrecision_dX = config.getDouble(ICSP::section, ICSP::integratedCrossSectionIntegralPrecision_dX);
    string approximationMethod = config.getValue(ICSP::section, ICSP::approximationMethod);
    int numberOfThreads = config.getInt(ICSP::section, ICSP::numberOfThreads);
  
    cout << ICSP::numberOfPointsIntegratedCrossSections << " = " << numberOfPointsIntegratedCrossSections << endl;
    cout << ICSP::propagatorIntegralPrecision << " = " << propagatorIntegralPrecision << endl;
    cout << ICSP::largeAngleScatteringContribution << " = " << largeAngleScatteringContribution << endl;
    cout << ICSP::crossSectionIntegralPrecision << " = " << crossSectionIntegralPrecision << endl;
    cout << ICSP::integratedCrossSectionIntegralPrecision_dXdY << " = " << integratedCrossSectionIntegralPrecision_dXdY << endl;
    cout << ICSP::integratedCrossSectionIntegralPrecision_dX << " = " << integratedCrossSectionIntegralPrecision_dX << endl;
    cout << ICSP::approximationMethod << " = " << approximationMethod << endl;
    cout << ICSP::numberOfThreads << " = " << numberOfThreads << endl;

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
    
    evaluateIsospinSymmetricIntegratedCrossSectionsWithZeroChemicalPotential(
        parameters, 
        precisionVacuum, 
        stringToMultiRootFindingMethod(methodVacuum), 
        upQuarkMassGuess, 
        strangeQuarkMassGuess, 
        precisionVacToFinTemp,
        stringToMultiRootFindingMethod(methodVacToFinTemp), 
        nearVacuumTemperature, 
        minimumTemp, 
        maximumTemp, 
        numberOfPointsFromVacToFinTemp, 
        numberOfPointsFromLowToHighTemp, 
        largeAngleScatteringContribution, 
        stringToIntegratedCrossSectionApproximationMethod(approximationMethod),
        propagatorIntegralPrecision,
        crossSectionIntegralPrecision,
        integratedCrossSectionIntegralPrecision_dXdY,
        integratedCrossSectionIntegralPrecision_dX,
        numberOfThreads
    );
}

bool IsospinSymmetricIntegratedCrossSectionsFiniteChemPot::validateTemperatureAndGridConsistency() const
{
    namespace VFTZCPP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    namespace LHT = SU3NJL3DCutoffFileParserKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters;
    namespace ICS = SU3NJL3DCutoffFileParserKeys::IntegratedCrossSectionsParameters;

    bool areTemperaturesEqual = (config.getValue(VFTZCPP::section, VFTZCPP::temperature)==config.getValue(LHT::section, LHT::minimumTemp));
    if ( !areTemperaturesEqual )
    {
        cout << invalidFileMessage + " Invalid values found in sections " + VFTZCPP::section + " and " + LHT::section + ": \n";
        cout << VFTZCPP::temperature << "==" << LHT::minimumTemp << " must be satisfied.\n";
    }

    bool areNumberOfPointsEqual = (config.getValue(LHT::section, LHT::numberOfPointsFromLowToHighTemp)==config.getValue(ICS::section, ICS::numberOfPointsIntegratedCrossSections));
    if ( !areNumberOfPointsEqual )
    {
        cout << invalidFileMessage + " Invalid values found in sections " + LHT::section + " and " + ICS::section + ": \n";
        cout << LHT::numberOfPointsFromLowToHighTemp << "==" << ICS::numberOfPointsIntegratedCrossSections << " must be satisfied.\n";
    }
    
    return areTemperaturesEqual && areNumberOfPointsEqual;
}

bool IsospinSymmetricIntegratedCrossSectionsFiniteChemPot::validateFile() const
{   
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
	const SU3NJL3DCutoffFileParser::Vacuum::Masses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    bool allRequiredSectionsPresent = true;
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters::section,
        SU3NJL3DCutoffFileParserKeys::FiniteTemperatureToFiniteChemicalPotentialParameters::section,
        SU3NJL3DCutoffFileParserKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters::section,
        SU3NJL3DCutoffFileParserKeys::IntegratedCrossSectionsParameters::section
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
    bool areLowToHighTemperatureAtZeroChemicalPotentialParametersValid = validateLowToHighTemperatureAtZeroChemicalPotentialParameters();
    bool areIntegratedCrossSectionsParametersValid = validateIntegratedCrossSectionsParameters();

    // Other validations
    bool areTemperatureAndGridValid = validateTemperatureAndGridConsistency();

    return vacuumValidations && 
           allRequiredSectionsPresent &&
           areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid &&
           areFiniteTemperatureToFiniteChemicalPotentialParametersValid &&
           areLowToHighTemperatureAtZeroChemicalPotentialParametersValid &&
           areIntegratedCrossSectionsParametersValid &&
           areTemperatureAndGridValid;
}

void IsospinSymmetricIntegratedCrossSectionsFiniteChemPot::evaluate() const
{   
    namespace MP = SU3NJL3DCutoffFileParserKeys::ModelParameters;
    namespace VMP = SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters;
    namespace VFTZCPP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    namespace FTFCPP = SU3NJL3DCutoffFileParserKeys::FiniteTemperatureToFiniteChemicalPotentialParameters;
    namespace LHT = SU3NJL3DCutoffFileParserKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters;
    namespace ICSP = SU3NJL3DCutoffFileParserKeys::IntegratedCrossSectionsParameters;

    // Model Parameters
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);

    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toString(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMassesParameters
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toString(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;
    
    // VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters
    cout << "\n" << VFTZCPP::section << ": " << endl;

    double nearVacuumTemperature = config.getDouble(VFTZCPP::section, VFTZCPP::nearVacuumTemperature);
    double temperature = config.getDouble(VFTZCPP::section, VFTZCPP::temperature);
    int numberOfPointsFromVacToFinTemp = config.getInt(VFTZCPP::section, VFTZCPP::numberOfPointsFromVacToFinTemp);
    double precisionVacToFinTemp = config.getDouble(VFTZCPP::section, VFTZCPP::precisionVacToFinTemp);
    string methodVacToFinTemp = config.getValue(VFTZCPP::section, VFTZCPP::methodVacToFinTemp);

    cout << VFTZCPP::nearVacuumTemperature << " = " << nearVacuumTemperature << endl;
    cout << VFTZCPP::temperature << " = " << temperature << endl;
    cout << VFTZCPP::numberOfPointsFromVacToFinTemp << " = " << numberOfPointsFromVacToFinTemp << endl;
    cout << VFTZCPP::precisionVacToFinTemp << " = " << precisionVacToFinTemp << endl;
    cout << VFTZCPP::methodVacToFinTemp << " = " << methodVacToFinTemp << endl;

    // FiniteTemperatureToFiniteChemicalPotentialParameters
    cout << "\n" << FTFCPP::section << ": " << endl;
    
    double chemPot = config.getDouble(FTFCPP::section, FTFCPP::chemPot);
    int numberOfPointsFromFinTempToFinChemPot = config.getInt(FTFCPP::section, FTFCPP::numberOfPointsFromFinTempToFinChemPot);
    double precisionFinTempToFinChemPot = config.getDouble(FTFCPP::section, FTFCPP::precisionFinTempToFinChemPot);
    string methodFinTempToFinChemPot = config.getValue(FTFCPP::section, FTFCPP::methodFinTempToFinChemPot);

    cout << FTFCPP::chemPot << " = " << chemPot << endl;
    cout << FTFCPP::numberOfPointsFromFinTempToFinChemPot << " = " << numberOfPointsFromFinTempToFinChemPot << endl;
    cout << FTFCPP::precisionFinTempToFinChemPot << " = " << precisionFinTempToFinChemPot << endl;
    cout << FTFCPP::methodFinTempToFinChemPot << " = " << methodFinTempToFinChemPot << endl;
    
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
    cout << "\n" << ICSP::section << ": " << endl;
    
    int numberOfPointsIntegratedCrossSections = config.getInt(ICSP::section, ICSP::numberOfPointsIntegratedCrossSections);
    double propagatorIntegralPrecision = config.getDouble(ICSP::section, ICSP::propagatorIntegralPrecision);
    bool largeAngleScatteringContribution = config.getBool(ICSP::section, ICSP::largeAngleScatteringContribution);
    double crossSectionIntegralPrecision = config.getDouble(ICSP::section, ICSP::crossSectionIntegralPrecision);
    double integratedCrossSectionIntegralPrecision_dXdY = config.getDouble(ICSP::section, ICSP::integratedCrossSectionIntegralPrecision_dXdY);
    double integratedCrossSectionIntegralPrecision_dX = config.getDouble(ICSP::section, ICSP::integratedCrossSectionIntegralPrecision_dX);
    string approximationMethod = config.getValue(ICSP::section, ICSP::approximationMethod);
    int numberOfThreads = config.getInt(ICSP::section, ICSP::numberOfThreads);
  
    cout << ICSP::numberOfPointsIntegratedCrossSections << " = " << numberOfPointsIntegratedCrossSections << endl;
    cout << ICSP::propagatorIntegralPrecision << " = " << propagatorIntegralPrecision << endl;
    cout << ICSP::largeAngleScatteringContribution << " = " << largeAngleScatteringContribution << endl;
    cout << ICSP::crossSectionIntegralPrecision << " = " << crossSectionIntegralPrecision << endl;
    cout << ICSP::integratedCrossSectionIntegralPrecision_dXdY << " = " << integratedCrossSectionIntegralPrecision_dXdY << endl;
    cout << ICSP::integratedCrossSectionIntegralPrecision_dX << " = " << integratedCrossSectionIntegralPrecision_dX << endl;
    cout << ICSP::approximationMethod << " = " << approximationMethod << endl;
    cout << ICSP::numberOfThreads << " = " << numberOfThreads << endl;

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
    
    evaluateIsospinSymmetricIntegratedCrossSectionsWithFixedChemicalPotential(
        parameters, 
        precisionVacuum, 
        stringToMultiRootFindingMethod(methodVacuum), 
        upQuarkMassGuess, 
        strangeQuarkMassGuess,
        nearVacuumTemperature,
        minimumTemp,
        numberOfPointsFromVacToFinTemp, 
        precisionVacToFinTemp, 
        stringToMultiRootFindingMethod(methodVacToFinTemp), 
        chemPot,
        numberOfPointsFromFinTempToFinChemPot,
        precisionFinTempToFinChemPot, 
        stringToMultiRootFindingMethod(methodFinTempToFinChemPot),
        maximumTemp, 
        numberOfPointsFromLowToHighTemp, 
        precisionLowToHighTemp, 
        stringToMultiRootFindingMethod(methodLowToHighTemp),
        largeAngleScatteringContribution, 
        stringToIntegratedCrossSectionApproximationMethod(approximationMethod),
        propagatorIntegralPrecision,
        crossSectionIntegralPrecision,
        integratedCrossSectionIntegralPrecision_dXdY,
        integratedCrossSectionIntegralPrecision_dX,
        numberOfThreads
    );
}

bool ThermoFixedChemPotTrajectory::validateFile() const
{   
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters
	const SU3NJL3DCutoffFileParser::Vacuum::Masses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffFileParserKeys::VacuumToFiniteChemicalPotentialParameters::section,
        SU3NJL3DCutoffFileParserKeys::ToTemperatureParameters::section,
        SU3NJL3DCutoffFileParserKeys::OutputFileParameters::section,
    };
    bool allRequiredSectionsPresent = checkRequiredSections(requiredSections);

    // Validate individual sections
    bool areVacuumToFiniteChemicalPotentialParametersValid = validateVacuumToFiniteChemicalPotentialParameters();
    bool areToTemperatureParametersValid = validateToTemperatureParameters();

    return vacuumValidations && 
           allRequiredSectionsPresent &&
           areVacuumToFiniteChemicalPotentialParametersValid &&
           areToTemperatureParametersValid;
}

void ThermoFixedChemPotTrajectory::evaluate() const
{   
    namespace MP = SU3NJL3DCutoffFileParserKeys::ModelParameters;
    namespace VMP = SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters;
    namespace VTFCPP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteChemicalPotentialParameters;
    namespace TTP = SU3NJL3DCutoffFileParserKeys::ToTemperatureParameters;
    namespace OFP = SU3NJL3DCutoffFileParserKeys::OutputFileParameters;

    // Model Parameters
    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);

    cout << "\n" << MP::section << ": " << endl;
    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toString(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMassesParameters
    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << "\n" << VMP::section << ": " << endl;
    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toString(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;
    
    // VacuumToFiniteChemicalPotentialParameters
    double chemicalPotential = config.getDouble(VTFCPP::section, VTFCPP::chemicalPotential);
    int numberOfPointsVacToChemPot = config.getInt(VTFCPP::section, VTFCPP::numberOfPointsVacToChemPot);
    double precisionVacToChemPot = config.getDouble(VTFCPP::section, VTFCPP::precisionVacToChemPot);
    string methodVacToChemPot = config.getValue(VTFCPP::section, VTFCPP::methodVacToChemPot);

    cout << "\n" << VTFCPP::section << ": " << endl;
    cout << VTFCPP::chemicalPotential << " = " << chemicalPotential << endl;
    cout << VTFCPP::numberOfPointsVacToChemPot << " = " << numberOfPointsVacToChemPot << endl;
    cout << VTFCPP::precisionVacToChemPot << " = " << precisionVacToChemPot << endl;
    cout << VTFCPP::methodVacToChemPot << " = " << methodVacToChemPot << endl;
    
    // ToTemperatureParameters
    double temperature = config.getDouble(TTP::section, TTP::temperature);
    int numberOfPointsUpToTemp = config.getInt(TTP::section, TTP::numberOfPointsUpToTemp);
	double precisionUpToTemp = config.getDouble(TTP::section, TTP::precisionUpToTemp);
	string methodUpToTemp = config.getValue(TTP::section, TTP::methodUpToTemp);

    cout << "\n" << TTP::section << ": " << endl;
    cout << TTP::temperature << " = " << temperature << endl;
    cout << TTP::numberOfPointsUpToTemp << " = " << numberOfPointsUpToTemp << endl;
    cout << TTP::precisionUpToTemp << " = " << precisionUpToTemp << endl;
    cout << TTP::methodUpToTemp << " = " << methodUpToTemp << endl;

    // OutputFileParameters    
    string customSuffix = config.getValue(OFP::section, OFP::customSuffix);
    
    cout << "\n" << OFP::section << ": " << endl;
    cout << OFP::customSuffix << " = " << customSuffix << endl;

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

    SU3NJL3DCutoffFixedChemPotTemp::computeThermoFixedChemPotTrajectory(
        parameters, 
        precisionVacuum,
        stringToMultiRootFindingMethod(methodVacuum),
        upQuarkMassGuess,
        downQuarkMassGuess,
        strangeQuarkMassGuess,
        chemicalPotential,
        numberOfPointsVacToChemPot,
        precisionVacToChemPot,
        stringToMultiRootFindingMethod(methodVacToChemPot),
        temperature,
        numberOfPointsUpToTemp,
        precisionUpToTemp,
        stringToMultiRootFindingMethod(methodUpToTemp),
        customSuffix
    );
}

bool ThermoFixedTemperatureTrajectory::validateFile() const
{   
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters
	const SU3NJL3DCutoffFileParser::Vacuum::Masses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffFileParserKeys::VacuumToTemperatureParameters::section,
        SU3NJL3DCutoffFileParserKeys::ToChemicalPotentialSymmetricParameters::section,
        SU3NJL3DCutoffFileParserKeys::OutputFileParameters::section,
    };
    bool allRequiredSectionsPresent = checkRequiredSections(requiredSections);

    // Validate individual sections
    bool areVacuumToTemperatureParametersValid = validateVacuumToTemperatureParameters();
    bool areToChemicalPotentialSymmetricParametersValid = validateToChemicalPotentialSymmetricParameters();

    return vacuumValidations && 
           allRequiredSectionsPresent &&
           areVacuumToTemperatureParametersValid &&
           areToChemicalPotentialSymmetricParametersValid;
}

void ThermoFixedTemperatureTrajectory::evaluate() const
{   
    namespace MP = SU3NJL3DCutoffFileParserKeys::ModelParameters;
    namespace VMP = SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters;
    namespace VTTP = SU3NJL3DCutoffFileParserKeys::VacuumToTemperatureParameters;
    namespace TCPSP = SU3NJL3DCutoffFileParserKeys::ToChemicalPotentialSymmetricParameters;
    namespace OFP = SU3NJL3DCutoffFileParserKeys::OutputFileParameters;

    // Model Parameters
    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);

    cout << "\n" << MP::section << ": " << endl;
    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toString(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMassesParameters
    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << "\n" << VMP::section << ": " << endl;
    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toString(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;

    // VacuumToTemperatureParameters
    double temperature = config.getDouble(VTTP::section, VTTP::temperature);
    int numberOfPointsVacToTemp = config.getInt(VTTP::section, VTTP::numberOfPointsVacToTemp);
    double precisionVacToTemp = config.getDouble(VTTP::section, VTTP::precisionVacToTemp);
    string methodVacToTemp = config.getValue(VTTP::section, VTTP::methodVacToTemp);

    cout << "\n" << VTTP::section << ": " << endl;
    cout << VTTP::temperature << " = " << temperature << endl;
    cout << VTTP::numberOfPointsVacToTemp << " = " << numberOfPointsVacToTemp << endl;
    cout << VTTP::precisionVacToTemp << " = " << precisionVacToTemp << endl;
    cout << VTTP::methodVacToTemp << " = " << methodVacToTemp << endl;
    
    // ToChemicalPotentialSymmetricParameters
    double chemicalPotential = config.getDouble(TCPSP::section, TCPSP::chemicalPotential);
    int numberOfPointsToChemPot = config.getInt(TCPSP::section, TCPSP::numberOfPointsToChemPot);
	double precisionToChemPot = config.getDouble(TCPSP::section, TCPSP::precisionToChemPot);
	string methodToChemPot = config.getValue(TCPSP::section, TCPSP::methodToChemPot);

    cout << "\n" << TCPSP::section << ": " << endl;
    cout << TCPSP::chemicalPotential << " = " << chemicalPotential << endl;
    cout << TCPSP::numberOfPointsToChemPot << " = " << numberOfPointsToChemPot << endl;
    cout << TCPSP::precisionToChemPot << " = " << precisionToChemPot << endl;
    cout << TCPSP::methodToChemPot << " = " << methodToChemPot << endl;

    // OutputFileParameters    
    string customSuffix = config.getValue(OFP::section, OFP::customSuffix);
    
    cout << "\n" << OFP::section << ": " << endl;
    cout << OFP::customSuffix << " = " << customSuffix << endl;

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

    SU3NJL3DCutoffFixedChemPotTemp::computeThermoFixedTemperatureTrajectory(
        parameters, 
        precisionVacuum,
        stringToMultiRootFindingMethod(methodVacuum),
        upQuarkMassGuess,
        downQuarkMassGuess,
        strangeQuarkMassGuess,
        temperature,
        numberOfPointsVacToTemp,
        precisionVacToTemp,
        stringToMultiRootFindingMethod(methodVacToTemp),
        chemicalPotential,
        numberOfPointsToChemPot,
        precisionToChemPot,
        stringToMultiRootFindingMethod(methodToChemPot),
        customSuffix
    );
}

bool InMediumMassesAndThermodynamics::validateFile() const
{   
    // Validate sections SU3NJL3DCutoffModelParameters, NJLDimensionfulCouplings and VacuumMassesParameters using previous developed logic
	const SU3NJL3DCutoffFileParser::Vacuum::Masses configVacuum(config);
    bool vacuumValidations = configVacuum.validateFile();

    // Check for missing sections
    vector<string> requiredSections = 
    {
        SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters::section
    };
    bool allRequiredSectionsPresent = checkRequiredSections(requiredSections);

    // Validate individual sections
    bool areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid = 
    validateVacuumToFiniteTemperatureAtZeroChemicalPotentialParameters();

    return vacuumValidations && 
           allRequiredSectionsPresent &&
           areVacuumToFiniteTemperatureAtZeroChemicalPotentialParametersValid;
}

void InMediumMassesAndThermodynamics::evaluate() const
{   
    namespace MP = SU3NJL3DCutoffFileParserKeys::ModelParameters;
    namespace VMP = SU3NJL3DCutoffFileParserKeys::VacuumMassesParameters;
    namespace VFTZCPP = SU3NJL3DCutoffFileParserKeys::VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters;
    namespace LHT = SU3NJL3DCutoffFileParserKeys::LowToHighTemperatureAtZeroChemicalPotentialParameters;
    namespace ICSP = SU3NJL3DCutoffFileParserKeys::IntegratedCrossSectionsParameters;

    // Model Parameters
    cout << "\n" << MP::section << ": " << endl;

    string parameterSetName = config.getValue(MP::section, MP::parameterSetName);
    string regularizationScheme = config.getValue(MP::section, MP::regularizationScheme);
    double cutoff = config.getDouble(MP::section, MP::cutoff);
    double upQuarkCurrentMass = config.getDouble(MP::section, MP::upQuarkCurrentMass);
    double downQuarkCurrentMass = config.getDouble(MP::section, MP::downQuarkCurrentMass);
    double strangeQuarkCurrentMass = config.getDouble(MP::section, MP::strangeQuarkCurrentMass);

    cout << MP::parameterSetName << " = " << parameterSetName << endl;
    cout << MP::regularizationScheme << " = " << toString(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << MP::cutoff << " = " << cutoff << endl;
    cout << MP::upQuarkCurrentMass << " = " << upQuarkCurrentMass << endl;
    cout << MP::downQuarkCurrentMass << " = " << downQuarkCurrentMass << endl;
    cout << MP::strangeQuarkCurrentMass << " = " << strangeQuarkCurrentMass << endl;

    // Dimensionful Couplings
    cout << "\n" << SU3NJL3DCutoffFileParserKeys::DimensionfulCouplings::section << ": " << endl;
    NJLDimensionfulCouplings couplings = extractDimensionfulCouplings();

    // VacuumMassesParameters
    cout << "\n" << VMP::section << ": " << endl;

    double precisionVacuum = config.getDouble(VMP::section, VMP::precisionVacuum);
    string methodVacuum = config.getValue(VMP::section, VMP::methodVacuum);
    double upQuarkMassGuess = config.getDouble(VMP::section, VMP::upQuarkMassGuess);
    double downQuarkMassGuess = config.getDouble(VMP::section, VMP::downQuarkMassGuess);
    double strangeQuarkMassGuess = config.getDouble(VMP::section, VMP::strangeQuarkMassGuess);

    cout << VMP::precisionVacuum << " = " << precisionVacuum << endl;
    cout << VMP::methodVacuum << " = " << toString(stringToMultiRootFindingMethod(methodVacuum)) << endl;
    cout << VMP::upQuarkMassGuess << " = " << upQuarkMassGuess << endl;
    cout << VMP::downQuarkMassGuess << " = " << downQuarkMassGuess << endl;
    cout << VMP::strangeQuarkMassGuess << " = " << strangeQuarkMassGuess << endl;
    
    // VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters
    cout << "\n" << VFTZCPP::section << ": " << endl;

    double nearVacuumTemperature = config.getDouble(VFTZCPP::section, VFTZCPP::nearVacuumTemperature);
    double temperature = config.getDouble(VFTZCPP::section, VFTZCPP::temperature);
    int numberOfPointsFromVacToFinTemp = config.getInt(VFTZCPP::section, VFTZCPP::numberOfPointsFromVacToFinTemp);
    double precisionVacToFinTemp = config.getDouble(VFTZCPP::section, VFTZCPP::precisionVacToFinTemp);
    string methodVacToFinTemp = config.getValue(VFTZCPP::section, VFTZCPP::methodVacToFinTemp);

    cout << VFTZCPP::nearVacuumTemperature << " = " << nearVacuumTemperature << endl;
    cout << VFTZCPP::temperature << " = " << temperature << endl;
    cout << VFTZCPP::numberOfPointsFromVacToFinTemp << " = " << numberOfPointsFromVacToFinTemp << endl;
    cout << VFTZCPP::precisionVacToFinTemp << " = " << precisionVacToFinTemp << endl;
    cout << VFTZCPP::methodVacToFinTemp << " = " << methodVacToFinTemp << endl;

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

    SU3NJL3DCutoffFixedChemPotTemp::evaluateInMediumMassesAndThermodynamics(
        parameters,                                    
        precisionVacuum,                                    
        stringToMultiRootFindingMethod(methodVacuum),                                  
        upQuarkMassGuess, 
        downQuarkMassGuess,
        strangeQuarkMassGuess,
        nearVacuumTemperature,
        temperature, 
        numberOfPointsFromVacToFinTemp,
        precisionVacToFinTemp,
        stringToMultiRootFindingMethod(methodVacToFinTemp)
    );
}

}

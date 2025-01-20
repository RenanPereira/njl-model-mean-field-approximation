#include "njl_model/n_fermion_line_integrals/KlevanskyB0Integral3DCutoffFileParser.h"
#include "njl_model/njl_regularization_schemes.h"


bool KlevanskyB0Integral3DCutoffFileParser::validateFileQuality() const 
{   
    // Check for missing sections
    vector<string> requiredSections = 
    {
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::section, 
    }; 
    // Validate sections
    if (!config.validateRequiredSections(requiredSections)) 
    {
        return false;
    }

    // Validate individual sections
    bool areVsK0ParametersValid = validateSectionVsK0Parameters();

    // The function return true only if all tests passed
    if ( areVsK0ParametersValid )
    { 
        return true; 
    }
    else
    { 
        return false; 
    }
}

bool KlevanskyB0Integral3DCutoffFileParser::validateSectionVsK0Parameters() const 
{   
    string section = KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::section;
    vector<string> requiredKeys = 
    {
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::numberOfPoints, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::k0LambdaRatioMin, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::k0LambdaRatioMax, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::regularizationScheme, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::temperature, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential1, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential2, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass1, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass2, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::threeMomentum, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::integralPrecision, 
    };

    // Check required keys
    if (!config.validateRequiredKeys(section, requiredKeys)) 
    {
        return false;
    }

    // Ensure numberOfPoints>0
    if (!config.validatePositiveInteger(
            section, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::numberOfPoints, 
            invalidFileMessage, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::numberOfPoints + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Validate regularizationScheme
    string regularizationScheme = 
        config.getValue(KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::section, 
                        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::regularizationScheme);
    if (!isValidNJL3DCutoffRegularizationScheme(regularizationScheme)) 
    {
        cout << invalidFileMessage << endl;
        return false;
    }

    // Ensure temperature>0
    if (!config.validateNonNegativeDouble(
            section, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::temperature, 
            invalidFileMessage, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::temperature + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Ensure effectiveChemicalPotential1>0
    if (!config.validateNonNegativeDouble(
            section, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential1, 
            invalidFileMessage, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential1 + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Ensure effectiveChemicalPotential2>0
    if (!config.validateNonNegativeDouble(
            section, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential2, 
            invalidFileMessage, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential2 + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Ensure effectiveMass1>0
    if (!config.validateNonNegativeDouble(
            section, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass1, 
            invalidFileMessage, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass1 + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Ensure effectiveMass2>0
    if (!config.validateNonNegativeDouble(
            section, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass2, 
            invalidFileMessage, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass2 + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Ensure threeMomentum>0
    if (!config.validateNonNegativeDouble(
            section, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::threeMomentum, 
            invalidFileMessage, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::threeMomentum + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Ensure integralPrecision>0
    if (!config.validateNonNegativeDouble(
            section, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::integralPrecision, 
            invalidFileMessage, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::integralPrecision + " > 0 must be satisfied.")) 
    {
        return false;
    }

    return true; // All validations passed
}

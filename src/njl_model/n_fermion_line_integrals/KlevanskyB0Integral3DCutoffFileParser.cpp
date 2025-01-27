#include "njl_model/n_fermion_line_integrals/KlevanskyB0Integral3DCutoffFileParser.h"
#include "njl_model/njl_regularization_schemes.h"

using namespace std;


bool KlevanskyB0Integral3DCutoffFileParser::validateFileQuality() const 
{   
    // Check for missing sections
    vector<string> vsK0section = 
    {
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::vsK0section, 
    }; 
    vector<string> vsKsection = 
    {
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::vsAbsKsection, 
    }; 
    if (!config.validateRequiredSections(vsK0section) && !config.validateRequiredSections(vsKsection)) 
    {
        return false;
    }

    // Get all sections for consequent validation
    vector<string> sections = config.getSectionOrder();
    
    bool areVsK0ParametersValid = true;
    for (int i = 0; i < int(sections.size()); i++)
    {
        if ( sections[i].find(KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::vsK0section)==0 )
        {
            //cout << "Validating section: " << sections[i] << ".\n";
            areVsK0ParametersValid = areVsK0ParametersValid && validateSectionVsK0Parameters(sections[i]);
        }
    }

    bool areVsAbsKParametersValid = true;
    for (int i = 0; i < int(sections.size()); i++)
    {
        if ( sections[i].find(KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::vsAbsKsection)==0 )
        {
            //cout << "Validating section: " << sections[i] << ".\n";
            areVsAbsKParametersValid = areVsAbsKParametersValid && validateSectionVsAbsKParameters(sections[i]);
        }
    }

    // The function return true only if all tests passed
    if ( areVsK0ParametersValid &&
         areVsAbsKParametersValid )
    { 
        return true; 
    }
    else
    { 
        return false; 
    }
    
    return true;
}

bool KlevanskyB0Integral3DCutoffFileParser::validateSectionVsK0Parameters(string sectionName) const 
{   
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
    if (!config.validateRequiredKeys(sectionName, requiredKeys)) 
    {
        return false;
    }

    // Ensure numberOfPoints>0
    if (!config.validatePositiveInteger(
            sectionName, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::numberOfPoints, 
            invalidFileMessage + " Invalid value found in section " + sectionName + ".", 
            KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::numberOfPoints + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Validate regularizationScheme
    string regularizationScheme = config.getValue(sectionName, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::regularizationScheme);
    if (!isValidNJL3DCutoffRegularizationScheme(regularizationScheme)) 
    {
        cout << invalidFileMessage + " Invalid value found in section " + sectionName + "." << endl;
        return false;
    }

    // Ensure all other parameters > 0
    vector<string> nonNegativeKeys = 
    {
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::temperature, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential1, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential2, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass1, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass2, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::threeMomentum, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::integralPrecision, 
    };
    for (int i = 0; i < int(nonNegativeKeys.size()); i++)
    {
        if (!config.validateNonNegativeDouble(
            sectionName, 
            nonNegativeKeys[i], 
            invalidFileMessage + " Invalid value found in section " + sectionName + ".", 
            nonNegativeKeys[i] + " > 0 must be satisfied.")) 
        {
            return false;
        }
    }

    return true; // All validations passed
}

bool KlevanskyB0Integral3DCutoffFileParser::validateSectionVsAbsKParameters(string sectionName) const 
{   
    vector<string> requiredKeys = 
    {
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::numberOfPoints, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::absKLambdaRatioMin, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::absKLambdaRatioMax, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::regularizationScheme, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::temperature, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveChemicalPotential1, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveChemicalPotential2, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveMass1, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveMass2, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::zeroMomentum, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::integralPrecision, 
    };

    // Check required keys
    if (!config.validateRequiredKeys(sectionName, requiredKeys)) 
    {
        return false;
    }

    // Ensure numberOfPoints>0
    if (!config.validatePositiveInteger(
            sectionName, 
            KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::numberOfPoints, 
            invalidFileMessage + " Invalid value found in section " + sectionName + ".", 
            KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::numberOfPoints + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Validate regularizationScheme
    string regularizationScheme = config.getValue(sectionName, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::regularizationScheme);
    if (!isValidNJL3DCutoffRegularizationScheme(regularizationScheme)) 
    {
        cout << invalidFileMessage + " Invalid value found in section " + sectionName + "." << endl;
        return false;
    }

    // Ensure all other parameters > 0
    vector<string> nonNegativeKeys = 
    {
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::temperature, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveChemicalPotential1, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveChemicalPotential2, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveMass1, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveMass2, 
        KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::integralPrecision, 
    };
    for (int i = 0; i < int(nonNegativeKeys.size()); i++)
    {
        if (!config.validateNonNegativeDouble(
            sectionName, 
            nonNegativeKeys[i], 
            invalidFileMessage + " Invalid value found in section " + sectionName + ".", 
            nonNegativeKeys[i] + " > 0 must be satisfied.")) 
        {
            return false;
        }
    }

    return true; // All validations passed
}

void KlevanskyB0Integral3DCutoffFileParser::evaluateB0VSK0OrAbsK() const
{   
    vector<map<string, string>> vsK0ParametersData = config.getSectionsData(KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::vsK0section);
    for (int i = 0; i < int(vsK0ParametersData.size()); ++i) 
	{
        cout << "\n" << KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::vsK0section << ":" << endl;
        cout << endl;

		const map<string, string>& section = vsK0ParametersData[i];

        int numberOfPoints = config.getInt(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::numberOfPoints);
        double k0LambdaRatioMin = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::k0LambdaRatioMin);
        double k0LambdaRatioMax = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::k0LambdaRatioMax);
        string regularizationScheme = config.getValue(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::regularizationScheme);
        double temperature = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::temperature);
        double effectiveChemicalPotential1 = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential1);
        double effectiveChemicalPotential2 = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveChemicalPotential2);
        double effectiveMass1 = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass1);
        double effectiveMass2 = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::effectiveMass2);
        double threeMomentum = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::threeMomentum);
        double integralPrecision = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsK0Parameters::integralPrecision);

        cout << "numberOfPoints = " << numberOfPoints << endl;
        cout << "k0LambdaRatioMin = " << k0LambdaRatioMin << endl;
        cout << "k0LambdaRatioMax = " << k0LambdaRatioMax << endl;
        cout << "regularizationScheme = " << regularizationScheme << endl;
        cout << "temperature = " << temperature << endl;
        cout << "effectiveChemicalPotential1 = " << effectiveChemicalPotential1 << endl;
        cout << "effectiveChemicalPotential2 = " << effectiveChemicalPotential2 << endl;
        cout << "effectiveMass1 = " << effectiveMass1 << endl;
        cout << "effectiveMass2 = " << effectiveMass2 << endl;
        cout << "threeMomentum = " << threeMomentum << endl;
        cout << "integralPrecision = " << integralPrecision << endl;
        cout << endl;

        // Call the function that performs the calculations

    }

    vector<map<string, string>> vsAbsKParametersData = config.getSectionsData(KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::vsAbsKsection);
    for (int i = 0; i < int(vsAbsKParametersData.size()); ++i) 
	{
        cout << "\n" << KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::vsAbsKsection << ":" << endl;
        cout << endl;

		const map<string, string>& section = vsAbsKParametersData[i];

        int numberOfPoints = config.getInt(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::numberOfPoints);
        double absKLambdaRatioMin = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::absKLambdaRatioMin);
        double absKLambdaRatioMax = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::absKLambdaRatioMax);
        string regularizationScheme = config.getValue(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::regularizationScheme);
        double temperature = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::temperature);
        double effectiveChemicalPotential1 = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveChemicalPotential1);
        double effectiveChemicalPotential2 = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveChemicalPotential2);
        double effectiveMass1 = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveMass1);
        double effectiveMass2 = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::effectiveMass2);
        double zeroMomentum = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::zeroMomentum);
        double integralPrecision = config.getDouble(section, KlevanskyB0Integral3DCutoffConfigKeys::VsAbsKParameters::integralPrecision);

        cout << "numberOfPoints = " << numberOfPoints << endl;
        cout << "absKLambdaRatioMin = " << absKLambdaRatioMin << endl;
        cout << "absKLambdaRatioMax = " << absKLambdaRatioMax << endl;
        cout << "regularizationScheme = " << regularizationScheme << endl;
        cout << "temperature = " << temperature << endl;
        cout << "effectiveChemicalPotential1 = " << effectiveChemicalPotential1 << endl;
        cout << "effectiveChemicalPotential2 = " << effectiveChemicalPotential2 << endl;
        cout << "effectiveMass1 = " << effectiveMass1 << endl;
        cout << "effectiveMass2 = " << effectiveMass2 << endl;
        cout << "zeroMomentum = " << zeroMomentum << endl;
        cout << "integralPrecision = " << integralPrecision << endl;
        cout << endl;

        // Call the function that performs the calculations

    }
}
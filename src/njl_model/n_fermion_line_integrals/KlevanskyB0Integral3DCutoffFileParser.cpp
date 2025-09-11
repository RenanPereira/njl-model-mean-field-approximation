#include "njl_model/n_fermion_line_integrals/KlevanskyB0Integral3DCutoffFileParser.h"
#include "njl_model/njl_regularization_schemes.h"
#include "njl_model/n_fermion_line_integrals/n_fermion_line_integrals_calculator.h"

using namespace std;


bool KlevanskyB0Integral3DCutoffFileParser::validateFileQualityEvaluateIntegral() const 
{   
    // Check for missing sections
    vector<string> vsK0section = 
    {
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::vsK0section, 
    }; 
    vector<string> vsKsection = 
    {
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::vsKsection, 
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
        if ( sections[i].find(KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::vsK0section)==0 )
        {
            //cout << "Validating section: " << sections[i] << ".\n";
            areVsK0ParametersValid = areVsK0ParametersValid && validateSectionVsK0Parameters(sections[i]);
        }
    }

    bool areVsKParametersValid = true;
    for (int i = 0; i < int(sections.size()); i++)
    {
        if ( sections[i].find(KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::vsKsection)==0 )
        {
            //cout << "Validating section: " << sections[i] << ".\n";
            areVsKParametersValid = areVsKParametersValid && validateSectionVsKParameters(sections[i]);
        }
    }

    // The function return true only if all tests passed
    if ( areVsK0ParametersValid &&
         areVsKParametersValid )
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
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::numberOfPoints, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::k0LambdaRatioMin, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::k0LambdaRatioMax, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::regularizationScheme, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::temperature, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveChemicalPotential1, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveChemicalPotential2, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::threeMomentumCutoff,
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveMass1, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveMass2, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::threeMomentum, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::integralPrecision, 
    };

    // Check required keys
    if (!config.validateRequiredKeys(sectionName, requiredKeys)) 
    {
        return false;
    }

    // Ensure numberOfPoints>0
    if (!config.validatePositiveInteger(
            sectionName, 
            KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::numberOfPoints, 
            invalidFileMessage + " Invalid value found in section " + sectionName + ".", 
            KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::numberOfPoints + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Validate regularizationScheme
    string regularizationScheme = config.getValue(sectionName, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::regularizationScheme);
    if (!isValidNJL3DCutoffRegularizationScheme(regularizationScheme)) 
    {
        cout << invalidFileMessage + " Invalid value found in section " + sectionName + "." << endl;
        return false;
    }

    // Ensure all other parameters > 0
    vector<string> nonNegativeKeys = 
    {
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::temperature, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveChemicalPotential1, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveChemicalPotential2, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::threeMomentumCutoff,
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveMass1, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveMass2, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::threeMomentum, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::integralPrecision, 
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

bool KlevanskyB0Integral3DCutoffFileParser::validateSectionVsKParameters(string sectionName) const 
{   
    vector<string> requiredKeys = 
    {
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::numberOfPoints, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::absKLambdaRatioMin, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::absKLambdaRatioMax, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::regularizationScheme, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::temperature, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveChemicalPotential1, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveChemicalPotential2, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::threeMomentumCutoff, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveMass1, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveMass2, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::zeroMomentum, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::integralPrecision, 
    };

    // Check required keys
    if (!config.validateRequiredKeys(sectionName, requiredKeys)) 
    {
        return false;
    }

    // Ensure numberOfPoints>0
    if (!config.validatePositiveInteger(
            sectionName, 
            KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::numberOfPoints, 
            invalidFileMessage + " Invalid value found in section " + sectionName + ".", 
            KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::numberOfPoints + " > 0 must be satisfied.")) 
    {
        return false;
    }

    // Validate regularizationScheme
    string regularizationScheme = config.getValue(sectionName, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::regularizationScheme);
    if (!isValidNJL3DCutoffRegularizationScheme(regularizationScheme)) 
    {
        cout << invalidFileMessage + " Invalid value found in section " + sectionName + "." << endl;
        return false;
    }

    // Ensure all other parameters > 0
    vector<string> nonNegativeKeys = 
    {   
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::absKLambdaRatioMin, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::temperature, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveChemicalPotential1, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveChemicalPotential2, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::threeMomentumCutoff, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveMass1, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveMass2, 
        KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::integralPrecision, 
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

void KlevanskyB0Integral3DCutoffFileParser::evaluateKlevanskyB0Integral3DCutoff() const
{   
    vector<map<string, string>> vsK0ParametersData = config.getSectionsData(KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::vsK0section);
    for (int i = 0; i < int(vsK0ParametersData.size()); ++i) 
	{
        cout << "\n" << KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::vsK0section << ":" << endl;
        cout << endl;

		const map<string, string>& section = vsK0ParametersData[i];

        int numberOfPoints = config.getInt(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::numberOfPoints);
        double k0LambdaRatioMin = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::k0LambdaRatioMin);
        double k0LambdaRatioMax = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::k0LambdaRatioMax);
        string regularizationScheme = config.getValue(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::regularizationScheme);
        double temperature = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::temperature);
        double effectiveChemicalPotential1 = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveChemicalPotential1);
        double effectiveChemicalPotential2 = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveChemicalPotential2);
        double threeMomentumCutoff = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::threeMomentumCutoff);
        double effectiveMass1 = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveMass1);
        double effectiveMass2 = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::effectiveMass2);
        double threeMomentum = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::threeMomentum);
        double integralPrecision = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::integralPrecision);

        cout << "numberOfPoints = " << numberOfPoints << endl;
        cout << "k0LambdaRatioMin = " << k0LambdaRatioMin << endl;
        cout << "k0LambdaRatioMax = " << k0LambdaRatioMax << endl;
        cout << "regularizationScheme = " << regularizationScheme << endl;
        cout << "temperature = " << temperature << endl;
        cout << "effectiveChemicalPotential1 = " << effectiveChemicalPotential1 << endl;
        cout << "effectiveChemicalPotential2 = " << effectiveChemicalPotential2 << endl;
        cout << "threeMomentumCutoff = " << threeMomentumCutoff << endl;
        cout << "effectiveMass1 = " << effectiveMass1 << endl;
        cout << "effectiveMass2 = " << effectiveMass2 << endl;
        cout << "threeMomentum = " << threeMomentum << endl;
        cout << "integralPrecision = " << integralPrecision << endl;
        cout << endl;

        // Call the function that performs the calculations
        evaluateKlevanskyB0Integral3DCutoffVsZeroMomentumToFile(
            numberOfPoints, 
            k0LambdaRatioMin, 
            k0LambdaRatioMax,
            stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
            temperature,
            effectiveChemicalPotential1,
            effectiveChemicalPotential2,
            threeMomentumCutoff,
            effectiveMass1,
            effectiveMass2,
            threeMomentum,
            integralPrecision);
    }

    vector<map<string, string>> vsKParametersData = config.getSectionsData(KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::vsKsection);
    for (int i = 0; i < int(vsKParametersData.size()); ++i) 
	{
        cout << "\n" << KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::vsKsection << ":" << endl;
        cout << endl;

		const map<string, string>& section = vsKParametersData[i];

        int numberOfPoints = config.getInt(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::numberOfPoints);
        double absKLambdaRatioMin = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::absKLambdaRatioMin);
        double absKLambdaRatioMax = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::absKLambdaRatioMax);
        string regularizationScheme = config.getValue(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::regularizationScheme);
        double temperature = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::temperature);
        double effectiveChemicalPotential1 = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveChemicalPotential1);
        double effectiveChemicalPotential2 = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveChemicalPotential2);
        double threeMomentumCutoff = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsK0Parameters::threeMomentumCutoff);
        double effectiveMass1 = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveMass1);
        double effectiveMass2 = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::effectiveMass2);
        double zeroMomentum = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::zeroMomentum);
        double integralPrecision = config.getDouble(section, KlevanskyB0Integral3DCutoffFileParserKeys::VsKParameters::integralPrecision);

        cout << "numberOfPoints = " << numberOfPoints << endl;
        cout << "absKLambdaRatioMin = " << absKLambdaRatioMin << endl;
        cout << "absKLambdaRatioMax = " << absKLambdaRatioMax << endl;
        cout << "regularizationScheme = " << regularizationScheme << endl;
        cout << "temperature = " << temperature << endl;
        cout << "effectiveChemicalPotential1 = " << effectiveChemicalPotential1 << endl;
        cout << "effectiveChemicalPotential2 = " << effectiveChemicalPotential2 << endl;
        cout << "threeMomentumCutoff = " << threeMomentumCutoff << endl;
        cout << "effectiveMass1 = " << effectiveMass1 << endl;
        cout << "effectiveMass2 = " << effectiveMass2 << endl;
        cout << "zeroMomentum = " << zeroMomentum << endl;
        cout << "integralPrecision = " << integralPrecision << endl;
        cout << endl;

        // Call the function that performs the calculations
        evaluateKlevanskyB0Integral3DCutoffVsThreeMomentumToFile(
            numberOfPoints, 
            absKLambdaRatioMin, 
            absKLambdaRatioMax,
            stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
            temperature,
            effectiveChemicalPotential1,
            effectiveChemicalPotential2,
            threeMomentumCutoff,
            effectiveMass1,
            effectiveMass2,
            zeroMomentum,
            integralPrecision);
    }
}

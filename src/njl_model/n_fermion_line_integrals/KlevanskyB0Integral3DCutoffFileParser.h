#ifndef KLEVANSKYB0INTEGRAL3DCUTOFFFILPARSER_H
#define KLEVANSKYB0INTEGRAL3DCUTOFFFILPARSER_H

#include "ini_file_parser/IniFileParser.h"
#include "njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff.h"


namespace KlevanskyB0Integral3DCutoffConfigKeys 
{
    namespace CalculationType
    {
        const std::string evaluateIntegral = "evaluateKlevanskyB0Integral3DCutoff";
    }

    namespace VsK0Parameters 
    {
        const std::string vsK0section = "VsK0Parameters";
        const std::string numberOfPoints = "numberOfPoints";
        const std::string k0LambdaRatioMin = "k0LambdaRatioMin";
        const std::string k0LambdaRatioMax = "k0LambdaRatioMax";
        const std::string regularizationScheme = "regularizationScheme";
        const std::string temperature = "temperature";
        const std::string effectiveChemicalPotential1 = "effectiveChemicalPotential1";
        const std::string effectiveChemicalPotential2 = "effectiveChemicalPotential2";
        const std::string threeMomentumCutoff = "threeMomentumCutoff";
        const std::string effectiveMass1 = "effectiveMass1";
        const std::string effectiveMass2 = "effectiveMass2";
        const std::string threeMomentum = "threeMomentum";
        const std::string integralPrecision = "integralPrecision";
    }

    namespace VsKParameters 
    {
        const std::string vsKsection = "VsKParameters";
        const std::string numberOfPoints = "numberOfPoints";
        const std::string absKLambdaRatioMin = "absKLambdaRatioMin";
        const std::string absKLambdaRatioMax = "absKLambdaRatioMax";
        const std::string regularizationScheme = "regularizationScheme";
        const std::string temperature = "temperature";
        const std::string effectiveChemicalPotential1 = "effectiveChemicalPotential1";
        const std::string effectiveChemicalPotential2 = "effectiveChemicalPotential2";
        const std::string threeMomentumCutoff = "threeMomentumCutoff";
        const std::string effectiveMass1 = "effectiveMass1";
        const std::string effectiveMass2 = "effectiveMass2";
        const std::string zeroMomentum = "zeroMomentum";
        const std::string integralPrecision = "integralPrecision";
    }
}


class KlevanskyB0Integral3DCutoffFileParser
{
public:
    const IniFileParser& config;    
    std::string invalidFileMessage = "Error: Invalid configuration found in the " + config.getFilename() + " file.";

public:
    KlevanskyB0Integral3DCutoffFileParser(const IniFileParser& p) : config(p) {};

    bool validateSectionVsK0Parameters(std::string ) const;
    bool validateSectionVsKParameters(std::string ) const;

    bool validateFileQualityEvaluateIntegral() const;
    
    void evaluateKlevanskyB0Integral3DCutoff() const;
};

#endif

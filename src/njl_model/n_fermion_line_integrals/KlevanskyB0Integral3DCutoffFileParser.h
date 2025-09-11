#ifndef KLEVANSKYB0INTEGRAL3DCUTOFFFILPARSER_H
#define KLEVANSKYB0INTEGRAL3DCUTOFFFILPARSER_H

#include "ini_file_parser/IniFileParser.h"
#include "njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff.h"


namespace KlevanskyB0Integral3DCutoffFileParserKeys 
{
    namespace VsK0Parameters 
    {
        inline const std::string vsK0section = "VsK0Parameters";
        inline const std::string numberOfPoints = "numberOfPoints";
        inline const std::string k0LambdaRatioMin = "k0LambdaRatioMin";
        inline const std::string k0LambdaRatioMax = "k0LambdaRatioMax";
        inline const std::string regularizationScheme = "regularizationScheme";
        inline const std::string temperature = "temperature";
        inline const std::string effectiveChemicalPotential1 = "effectiveChemicalPotential1";
        inline const std::string effectiveChemicalPotential2 = "effectiveChemicalPotential2";
        inline const std::string threeMomentumCutoff = "threeMomentumCutoff";
        inline const std::string effectiveMass1 = "effectiveMass1";
        inline const std::string effectiveMass2 = "effectiveMass2";
        inline const std::string threeMomentum = "threeMomentum";
        inline const std::string integralPrecision = "integralPrecision";
    }

    namespace VsKParameters 
    {
        inline const std::string vsKsection = "VsKParameters";
        inline const std::string numberOfPoints = "numberOfPoints";
        inline const std::string absKLambdaRatioMin = "absKLambdaRatioMin";
        inline const std::string absKLambdaRatioMax = "absKLambdaRatioMax";
        inline const std::string regularizationScheme = "regularizationScheme";
        inline const std::string temperature = "temperature";
        inline const std::string effectiveChemicalPotential1 = "effectiveChemicalPotential1";
        inline const std::string effectiveChemicalPotential2 = "effectiveChemicalPotential2";
        inline const std::string threeMomentumCutoff = "threeMomentumCutoff";
        inline const std::string effectiveMass1 = "effectiveMass1";
        inline const std::string effectiveMass2 = "effectiveMass2";
        inline const std::string zeroMomentum = "zeroMomentum";
        inline const std::string integralPrecision = "integralPrecision";
    }
}


class KlevanskyB0Integral3DCutoffFileParser
{
    public:
        inline static const std::string klevanskyB0Integral3DCutoff = "KlevanskyB0Integral3DCutoff";

    public:
        const IniFileParser& config;    
        std::string invalidFileMessage;

    public:
        KlevanskyB0Integral3DCutoffFileParser(const IniFileParser& p)
            : config(p),
            invalidFileMessage("Error: Invalid configuration found in the " + p.getFilename() + " file.")
        {}

        bool validateSectionVsK0Parameters(std::string ) const;
        bool validateSectionVsKParameters(std::string ) const;
        bool validateFileQualityEvaluateIntegral() const;
        void evaluateKlevanskyB0Integral3DCutoff() const;
};

#endif

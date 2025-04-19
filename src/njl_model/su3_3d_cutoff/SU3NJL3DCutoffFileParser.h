#ifndef SU3NJL3DCUTOFFFILEPARSER_H
#define SU3NJL3DCUTOFFFILEPARSER_H

#include "ini_file_parser/IniFileParser.h"
#include "njl_model/NJLDimensionfulCouplings.h"


namespace SU3NJL3DCutoffConfigKeys 
{
    namespace CalculationType
    {
        const std::string evaluateVacuumMasses = "evaluateSU3NJL3DCutoffVacuumMasses";
        const std::string evaluateFirstOrderLine = "evaluateSU3NJL3DCutoffFirstOrderLine";
    }

    namespace ModelParameters 
    {
        const std::string section = "SU3NJL3DCutoffModelParameters";
        const std::string parameterSetName = "parameterSetName";
        const std::string regularizationScheme = "regularizationScheme";
        const std::string cutoff = "cutoff_GeV";
        const std::string upQuarkCurrentMass = "upQuarkCurrentMass_GeV";
        const std::string downQuarkCurrentMass = "downQuarkCurrentMass_GeV";
        const std::string strangeQuarkCurrentMass = "strangeQuarkCurrentMass_GeV";
    }

    namespace DimensionfulCouplings 
    {
        const std::string section = "NJLDimensionfulCouplings";
        const std::string lagrangianInteractions = "lagrangianInteractions";
    }

    namespace VacuumMassesParameters 
    {
        const std::string section = "VacuumMassesParameters";
        const std::string precisionVacuum = "precisionVacuum";
        const std::string methodVacuum = "methodVacuum";
        const std::string upQuarkMassGuess = "upQuarkMassGuess_GeV";
        const std::string downQuarkMassGuess = "downQuarkMassGuess_GeV";
        const std::string strangeQuarkMassGuess = "strangeQuarkMassGuess_GeV";
    }

    namespace VacuumToFiniteBaryonDensityParameters
    {
        const std::string section = "VacuumToFiniteBaryonDensityParameters";
        const std::string minimumBaryonDensity = "minimumBaryonDensity_fmMinus3";
        const std::string maximumBaryonDensity = "maximumBaryonDensity_fmMinus3";
        const std::string numberOfPoints = "numberOfPoints";
        const std::string precisionZeroTempSol = "precisionZeroTemperatureSolution";
        const std::string methodZeroTempSol= "methodZeroTemperatureSolution";
    }

    namespace FirstOrderLineParameters
    {
        const std::string section = "FirstOrderLineParameters";
        const std::string precisionTransitionPointSol = "precisionTransitionPointSolution";
        const std::string methodTransitionPointSol = "methodTransitionPointSolution";
        const std::string deltaT = "deltaTemperature_GeV";
        const std::string massDifferenceCEP = "massDifferenceCEP_GeV";
    }
}


class SU3NJL3DCutoffFileParser
{
public:
    const IniFileParser& config;    
    std::string invalidFileMessage = "Error: Invalid configuration found in the " + config.getFilename() + " file.";

public:
    SU3NJL3DCutoffFileParser(const IniFileParser& p) : config(p) {}
    
    // Validations
    static NJLDimensionfulCouplings extractDimensionfulCouplings(const IniFileParser& );
    
    bool validateModelParameters() const;
    bool validateDimensionfulCouplings() const;
    bool validateVacuumMassesParameters() const;

    bool validateVacuumToFiniteBaryonDensityParameters() const;
    bool validateFirstOrderLineParameters() const;

    bool validateFileQualityEvaluateVacuumMasses() const;
    bool validateFileQualityEvaluateFirstOrderLine() const;

    //Calculations
    void evaluateVacuumMasses() const;
    void evaluateFirstOrderLine() const;
};

#endif
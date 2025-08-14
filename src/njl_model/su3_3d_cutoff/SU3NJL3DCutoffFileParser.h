#ifndef SU3NJL3DCUTOFFFILEPARSER_H
#define SU3NJL3DCUTOFFFILEPARSER_H

#include "ini_file_parser/IniFileParser.h"
#include "njl_model/NJLDimensionfulCouplings.h"


namespace SU3NJL3DCutoffConfigKeys 
{
    namespace CalculationType
    {
        const std::string Vacuum_evaluateVacuumMasses = "SU3NJL3DCutoffVacuum_evaluateVacuumMasses";
        const std::string FixedTempRhoBEqualChemPot_evaluateFirstOrderLine = "SU3NJL3DCutoffFixedTemperatureRhoBEqualChemicalPotential_evaluateFirstOrderLine";
        const std::string FixedChemPotTemp_evaluateCrossSectionsEqualLightMasses = "SU3NJL3DCutoffFixedChemicalPotentialTemperature_evaluateCrossSectionsEqualLightMasses";
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

    namespace VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters 
    {
        const std::string section = "VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters";
        const std::string temperature = "temperature_GeV";
        const std::string numberOfPointsFromVacToFinTemp = "numberOfPointsFromVacuumToFineteTemperature";
        const std::string precisionVacToFinTemp = "precisionVacuumToFiniteTemperature";
        const std::string methodVacToFinTemp = "methodVacuumToFiniteTemperature";
    }

    namespace FiniteTemperatureToFiniteChemicalPotentialParameters
    {
        const std::string section = "FiniteTemperatureToFiniteChemicalPotentialParameters";
        const std::string chemPot = "quarkChemicalPotential_GeV";
        const std::string numberOfPointsFromFinTempToFinChemPot = "numberOfPointsFromFiniteTemperatureToFiniteChemicalPotential";
        const std::string precisionFinTempToFinChemPot = "precisionFiniteTemperatureToFiniteChemicalPotential";
        const std::string methodFinTempToFinChemPot = "methodFiniteTemperatureToFiniteChemicalPotential";
    }
    
    namespace CrossSectionsParameters
    {
        const std::string section = "CrossSectionsParameters";
        const std::string propagatorIntegralPrecision = "propagatorIntegralPrecision";
        const std::string largeAngleScatteringContribution = "largeAngleScatteringContribution";
        const std::string precisionCrossSections = "precisionCrossSections";
        const std::string numberOfPointsCrossSections = "numberOfPointsCrossSections";
    }
}


class SU3NJL3DCutoffVacuumFileParser
{
public:
    const IniFileParser& config;    
    std::string invalidFileMessage = "Error: Invalid configuration found in the " + config.getFilename() + " file.";

public:
    SU3NJL3DCutoffVacuumFileParser(const IniFileParser& p) : config(p) {}
    
    // Validations
    static NJLDimensionfulCouplings extractDimensionfulCouplings(const IniFileParser& );
    
    bool validateModelParameters() const;
    bool validateDimensionfulCouplings() const;
    bool validateVacuumMassesParameters() const;
    bool validateFileQualityEvaluateVacuumMasses() const;

    void evaluateVacuumMasses() const;
};


class SU3NJL3DCutoffFixedTempRhoBEqualChemPotFileParser : public SU3NJL3DCutoffVacuumFileParser
{
public:
    SU3NJL3DCutoffFixedTempRhoBEqualChemPotFileParser(const IniFileParser& p) : SU3NJL3DCutoffVacuumFileParser(p) {}

    bool validateVacuumToFiniteBaryonDensityParameters() const;
    bool validateFirstOrderLineParameters() const;
    bool validateFileQualityEvaluateFirstOrderLine() const;

    void evaluateFirstOrderLine() const;
};


class SU3NJL3DCutoffFixedChemPotTempCrossSectionsFileParser : public SU3NJL3DCutoffVacuumFileParser
{
public:
    SU3NJL3DCutoffFixedChemPotTempCrossSectionsFileParser(const IniFileParser& p) : SU3NJL3DCutoffVacuumFileParser(p) {}

    bool validateVacuumToFiniteTemperatureAtZeroChemicalPotentialParameters() const;
    bool validateFiniteTemperatureToFiniteChemicalPotentialParameters() const;
    bool validateCrossSectionsParameters() const;
    bool validateFileQualityEvaluateCrossSectionsEqualLightMasses() const;

    void evaluateCrossSectionsEqualLightMasses() const;
};

#endif
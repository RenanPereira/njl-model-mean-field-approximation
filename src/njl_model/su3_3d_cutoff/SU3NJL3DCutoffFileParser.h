#ifndef SU3NJL3DCUTOFFFILEPARSER_H
#define SU3NJL3DCUTOFFFILEPARSER_H

#include "ini_file_parser/IniFileParser.h"
#include "njl_model/NJLDimensionfulCouplings.h"

namespace SU3NJL3DCutoffFileParserKeys 
{
    namespace ModelParameters 
    {
        inline const std::string section = "SU3NJL3DCutoffModelParameters";
        inline const std::string parameterSetName = "parameterSetName";
        inline const std::string regularizationScheme = "regularizationScheme";
        inline const std::string cutoff = "cutoff_GeV";
        inline const std::string upQuarkCurrentMass = "upQuarkCurrentMass_GeV";
        inline const std::string downQuarkCurrentMass = "downQuarkCurrentMass_GeV";
        inline const std::string strangeQuarkCurrentMass = "strangeQuarkCurrentMass_GeV";
    }

    namespace DimensionfulCouplings 
    {
        inline const std::string section = "NJLDimensionfulCouplings";
        inline const std::string lagrangianInteractions = "lagrangianInteractions";
    }

    namespace VacuumMassesParameters 
    {
        inline const std::string section = "VacuumMassesParameters";
        inline const std::string precisionVacuum = "precisionVacuum";
        inline const std::string methodVacuum = "methodVacuum";
        inline const std::string upQuarkMassGuess = "upQuarkMassGuess_GeV";
        inline const std::string downQuarkMassGuess = "downQuarkMassGuess_GeV";
        inline const std::string strangeQuarkMassGuess = "strangeQuarkMassGuess_GeV";
    }

    namespace VacuumToFiniteBaryonDensityParameters
    {
        inline const std::string section = "VacuumToFiniteBaryonDensityParameters";
        inline const std::string minimumBaryonDensity = "minimumBaryonDensity_fmMinus3";
        inline const std::string maximumBaryonDensity = "maximumBaryonDensity_fmMinus3";
        inline const std::string numberOfPoints = "numberOfPoints";
        inline const std::string precisionZeroTempSol = "precisionZeroTemperatureSolution";
        inline const std::string methodZeroTempSol= "methodZeroTemperatureSolution";
    }

    namespace FirstOrderLineParameters
    {
        inline const std::string section = "FirstOrderLineParameters";
        inline const std::string precisionTransitionPointSol = "precisionTransitionPointSolution";
        inline const std::string methodTransitionPointSol = "methodTransitionPointSolution";
        inline const std::string deltaT = "deltaTemperature_GeV";
        inline const std::string massDifferenceCEP = "massDifferenceCEP_GeV";
    }

    namespace VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters 
    {
        inline const std::string section = "VacuumToFiniteTemperatureAtZeroChemicalPotentialParameters";
        inline const std::string nearVacuumTemperature = "nearVacuumTemperature_GeV";
        inline const std::string temperature = "temperature_GeV";
        inline const std::string numberOfPointsFromVacToFinTemp = "numberOfPointsFromVacuumToFiniteTemperature";
        inline const std::string precisionVacToFinTemp = "precisionVacuumToFiniteTemperature";
        inline const std::string methodVacToFinTemp = "methodVacuumToFiniteTemperature";
    }

    namespace LowToHighTemperatureAtZeroChemicalPotentialParameters 
    {
        inline const std::string section = "LowToHighTemperatureAtZeroChemicalPotentialParameters";
        inline const std::string minimumTemp = "minimumTemperature_GeV";
        inline const std::string maximumTemp = "maximumTemperature_GeV";
        inline const std::string numberOfPointsFromLowToHighTemp = "numberOfPointsFromLowToHighTemperature";
        inline const std::string precisionLowToHighTemp = "precisionLowToHighTemperature";
        inline const std::string methodLowToHighTemp = "methodLowToHighTemperature";
    }    

    namespace FiniteTemperatureToFiniteChemicalPotentialParameters
    {
        inline const std::string section = "FiniteTemperatureToFiniteChemicalPotentialParameters";
        inline const std::string chemPot = "quarkChemicalPotential_GeV";
        inline const std::string numberOfPointsFromFinTempToFinChemPot = "numberOfPointsFromFiniteTemperatureToFiniteChemicalPotential";
        inline const std::string precisionFinTempToFinChemPot = "precisionFiniteTemperatureToFiniteChemicalPotential";
        inline const std::string methodFinTempToFinChemPot = "methodFiniteTemperatureToFiniteChemicalPotential";
    }
    
    namespace CrossSectionsParameters
    {
        inline const std::string section = "CrossSectionsParameters";
        inline const std::string propagatorIntegralPrecision = "propagatorIntegralPrecision";
        inline const std::string largeAngleScatteringContribution = "largeAngleScatteringContribution";
        inline const std::string precisionCrossSections = "precisionCrossSections";
        inline const std::string numberOfPointsCrossSections = "numberOfPointsCrossSections";
        inline const std::string numberOfThreads = "numberOfThreads";
    }

    namespace IntegratedCrossSectionsParameters
    {
        inline const std::string section = "IntegratedCrossSectionsParameters";
        inline const std::string numberOfPointsIntegratedCrossSections = "numberOfPointsIntegratedCrossSections";        
        inline const std::string propagatorIntegralPrecision = "propagatorIntegralPrecision";
        inline const std::string largeAngleScatteringContribution = "largeAngleScatteringContribution";
        inline const std::string crossSectionIntegralPrecision = "crossSectionIntegralPrecision";
        inline const std::string integratedCrossSectionIntegralPrecision_dXdY = "integratedCrossSectionIntegralPrecision_dXdY";
        inline const std::string integratedCrossSectionIntegralPrecision_dX = "integratedCrossSectionIntegralPrecision_dX";
        inline const std::string approximationMethod = "approximationMethod";
        inline const std::string numberOfThreads = "numberOfThreads";
    }
}

namespace SU3NJL3DCutoffFileParser
{   
    inline const std::string model = "SU3NJL3DCutoff";

    class Common
    {
        public:
            const IniFileParser& config;    
            std::string invalidFileMessage;

            Common(const IniFileParser& p)
                : config(p),
                invalidFileMessage("Error: Invalid configuration found in the " + p.getFilename() + " file.") {}

            virtual ~Common() {}

            // methods that need to be defined in children classes
            virtual bool validateFile() const = 0;
            virtual void evaluate() const = 0;
            void run() const
            {
                if (validateFile())
                {
                    evaluate();
                }
                else
                {
                    printQualityCheckFailedMessage();
                }
            }

            NJLDimensionfulCouplings extractDimensionfulCouplings() const;

            bool validateModelParameters() const;
            bool validateDimensionfulCouplings() const;
            bool validateVacuumMassesParameters() const;
            bool validateVacuumToFiniteBaryonDensityParameters() const;
            bool validateFirstOrderLineParameters() const;
            bool validateVacuumToFiniteTemperatureAtZeroChemicalPotentialParameters() const;
            bool validateFiniteTemperatureToFiniteChemicalPotentialParameters() const;
            bool validateCrossSectionsParameters() const;
            bool validateLowToHighTemperatureAtZeroChemicalPotentialParameters() const;
            bool validateIntegratedCrossSectionsParameters() const;

            bool checkRequiredSections(const std::vector<std::string> ) const;
            void printQualityCheckFailedMessage() const;
    };

    namespace Vacuum
    {   
        inline const std::string type = model + "Vacuum";

        class VacuumMasses : public Common
        {   
            public:
                inline static const std::string calculationType = type + "VacuumMasses";
            
            public:
                VacuumMasses(const IniFileParser& p) : Common(p) {}

                bool validateFile() const;
                void evaluate() const;
        };
    }

    namespace FixedTempRhoBEqualChemPot
    {
        inline const std::string type = model + "FixedTemperatureRhoBEqualChemicalPotential";

        class FirstOrderLine : public Common
        {   
            public:
                inline static const std::string calculationType = type + "FirstOrderLine";

            public:
                FirstOrderLine(const IniFileParser& p) : Common(p) {}

                bool validateFile() const;
                void evaluate() const;
        };
    }

    namespace FixedChemPotTemp
    {   
        inline const std::string type = model + "FixedChemicalPotentialTemperature";

        class IsospinSymmetricCrossSections : public Common
        {
            public:
                inline static const std::string calculationType = type + "IsospinSymmetricCrossSections";

            public:
                IsospinSymmetricCrossSections(const IniFileParser& p) : Common(p) {}
                
                bool validateFile() const;
                void evaluate() const;
        };

        class IsospinSymmetricIntegratedCrossSectionsZeroChemPot : public Common
        {   
            public:
                inline static const std::string calculationType = type + "IsospinSymmetricIntegratedCrossSectionsZeroChemicalPotential";

            public:
                IsospinSymmetricIntegratedCrossSectionsZeroChemPot(const IniFileParser& p) : Common(p) {}
                
                bool validateTemperatureAndGridConsistency() const;
                bool validateFile() const;
                void evaluate() const;
        };

        class InMediumMassesAndThermodynamics : public Common
        {   
            public:
                inline static const std::string calculationType = type + "InMediumMassesAndThermodynamics";

            public:
                InMediumMassesAndThermodynamics(const IniFileParser& p) : Common(p) {}
                
                bool validateFile() const;
                void evaluate() const;
        };
    }
}

#endif
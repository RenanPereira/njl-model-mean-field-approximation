#ifndef SU3NJL3DCUTOFFFILEPARSER_H
#define SU3NJL3DCUTOFFFILEPARSER_H

#include "ini_file_parser/IniFileParser.h"
#include "njl_model/NJLDimensionfulCouplings.h"

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
            bool validateVacuumToFiniteChemicalPotentialParameters() const;
            bool validateVacuumToTemperatureParameters() const;
            bool validateToTemperatureParameters() const;
            bool validateToChemicalPotentialSymmetricParameters() const;
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

        class Masses : public Common
        {   
            public:
                inline static const std::string calculationType = type + "Masses";
            
            public:
                Masses(const IniFileParser& p) : Common(p) {}

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
    
        class IsospinSymmetricIntegratedCrossSectionsFiniteChemPot : public Common
        {   
            public:
                inline static const std::string calculationType = type + "IsospinSymmetricIntegratedCrossSectionsFiniteChemicalPotential";

            public:
                IsospinSymmetricIntegratedCrossSectionsFiniteChemPot(const IniFileParser& p) : Common(p) {}
                
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

        class ThermoFixedChemPotTrajectory : public Common
        {   
            public:
                inline static const std::string calculationType = type + "ThermodynamicsFixedChemicalPotentialTrajectory";

            public:
                ThermoFixedChemPotTrajectory(const IniFileParser& p) : Common(p) {}
                
                bool validateFile() const;
                void evaluate() const;
        };

        class ThermoFixedTemperatureTrajectory : public Common
        {   
            public:
                inline static const std::string calculationType = type + "ThermodynamicsFixedTemperatureTrajectory";

            public:
                ThermoFixedTemperatureTrajectory(const IniFileParser& p) : Common(p) {}
                
                bool validateFile() const;
                void evaluate() const;
        };

    }
}

#endif
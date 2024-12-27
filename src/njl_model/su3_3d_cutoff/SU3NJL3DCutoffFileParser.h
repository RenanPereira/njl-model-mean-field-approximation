#ifndef SU3NJL3DCUTOFFFILEPARSER_H
#define SU3NJL3DCUTOFFFILEPARSER_H

#include "ini_file_parser/IniFileParser.h"
#include "njl_model/NJLDimensionfulCouplings.h"


namespace SU3NJL3DCutoffConfigKeys 
{
    namespace ModelParameters 
    {
        const std::string section = "SU3NJL3DCutoffModelParameters";
        const std::string parameterSetName = "parameterSetName";
        const std::string regularizationScheme = "regularizationScheme";
        const std::string cutoffInGeV = "cutoffInGeV";
        const std::string upQuarkCurrentMassInGeV = "upQuarkCurrentMassInGeV";
        const std::string downQuarkCurrentMassInGeV = "downQuarkCurrentMassInGeV";
        const std::string strangeQuarkCurrentMassInGeV = "strangeQuarkCurrentMassInGeV";
    }

    namespace DimensionfulCouplings 
    {
        const std::string section = "NJLDimensionfulCouplings";
        const std::string lagrangianInteractions = "lagrangianInteractions";
    }

    namespace GapEquationsVacuumParameters 
    {
        const std::string section = "SU3NJL3DCutoffGapEquationsVacuumParameters";
        const std::string gapPrecision = "gapPrecision";
        const std::string rootFindingMethod = "rootFindingMethod";
        const std::string upQuarkMassGuess = "upQuarkMassGuess";
        const std::string downQuarkMassGuess = "downQuarkMassGuess";
        const std::string strangeQuarkMassGuess = "strangeQuarkMassGuess";
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
    bool validateFileQuality() const;
    bool validateModelParameters() const;
    bool validateDimensionfulCouplings() const;
    bool validateGapEquationsVacuumParameters() const;

    //Calculations
    void evaluateVacuumMasses() const;
};

#endif
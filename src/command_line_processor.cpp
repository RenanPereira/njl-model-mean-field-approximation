#include "command_line_processor.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "njl_model/n_fermion_line_integrals/KlevanskyB0Integral3DCutoffFileParser.h"


int commandLineArgsProcessor(int argc, char* argv[])
{	
	std::string use_file_command = "use-config-file";

	// Handle command line input
    if (argc > 1) 
	{
        std::string command = argv[1];
        // Check for the commands
        if (command == use_file_command) 
		{
            if (argc == 3) 
			{	
				// Get name of the file from the third argument
                std::string configFileName = argv[2];

				// Open configuration file and parse it
				std::cout << "\nFeeding IniFileParser with file " << configFileName << "..." << std::endl;
				IniFileParser configFile(configFileName);
				selectPathBasedOnFileDetails(configFile);
            } 
			else 
			{
                std::cerr << "Error: No file provided after the command: " << command << std::endl;
                return 1;
            }
        } 
		else 
		{
            std::cerr << "Unknown command: " << command << std::endl;
            return 1;
        }
    }

	return 0;
}

void selectPathBasedOnFileDetails(const IniFileParser& configFile)
{	
    // Get the file type
    std::string type = configFile.getValue("FileDetails", "type");
	
	std::cout << "\nFileDetails:" << std::endl;
	std::cout << "type = " << type << std::endl;

	// Check if file is written correctly and then make calculation
	if(type==SU3NJL3DCutoffFileParser::Vacuum::Masses::calculationType)
	{	
		const SU3NJL3DCutoffFileParser::Vacuum::Masses config(configFile);
		config.run();
	}
	else if(type==KlevanskyB0Integral3DCutoffFileParser::calculationType)
	{	
		const KlevanskyB0Integral3DCutoffFileParser config(configFile);
		config.run();
	}
	else if (type==SU3NJL3DCutoffFileParser::FixedTempRhoBEqualChemPot::FirstOrderLine::calculationType)
	{
		const SU3NJL3DCutoffFileParser::FixedTempRhoBEqualChemPot::FirstOrderLine config(configFile);
		config.run();
	}
	else if (type==SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricCrossSections::calculationType)
	{
		const SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricCrossSections config(configFile);
		config.run();
	}
	else if (type==SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricIntegratedCrossSectionsZeroChemPot::calculationType)
	{	
		const SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricIntegratedCrossSectionsZeroChemPot config(configFile);
		config.run();
	}
	else if (type==SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricIntegratedCrossSectionsFiniteChemPot::calculationType)
	{	
		const SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricIntegratedCrossSectionsFiniteChemPot config(configFile);
		config.run();
	}
	else if (type==SU3NJL3DCutoffFileParser::FixedChemPotTemp::InMediumMassesAndThermodynamics::calculationType)
	{	
		const SU3NJL3DCutoffFileParser::FixedChemPotTemp::InMediumMassesAndThermodynamics config(configFile);
		config.run();
	}
	else if (type==SU3NJL3DCutoffFileParser::FixedChemPotTemp::ThermoFixedChemPotTrajectory::calculationType)
	{	
		const SU3NJL3DCutoffFileParser::FixedChemPotTemp::ThermoFixedChemPotTrajectory config(configFile);
		config.run();
	}
	else if (type==SU3NJL3DCutoffFileParser::FixedChemPotTemp::ThermoFixedTemperatureTrajectory::calculationType)
	{	
		const SU3NJL3DCutoffFileParser::FixedChemPotTemp::ThermoFixedTemperatureTrajectory config(configFile);
		config.run();
	}
	else
	{
		std::cout << "The file " << configFile.getFilename() << " does not match any known configuration! Check the FileDetails.\n";
	}
}

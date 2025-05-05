#include "command_line_processor.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "njl_model/n_fermion_line_integrals/KlevanskyB0Integral3DCutoffFileParser.h"

using namespace std;


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
				cout << "\nFeeding IniFileParser with file " << configFileName << "..." << endl;
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
    string fileTypeStr = configFile.getValue("FileDetails", "type");
	
	cout << "\nFileDetails:" << endl;
	cout << "type = " << fileTypeStr << endl;

	if(fileTypeStr==SU3NJL3DCutoffConfigKeys::CalculationType::Vacuum_evaluateVacuumMasses)
	{	
		//Check if file is written correctly
		const SU3NJL3DCutoffVacuumFileParser config(configFile);
		if(config.validateFileQualityEvaluateVacuumMasses())
		{	
			config.evaluateVacuumMasses();
		}
		else
		{
			cout << "The quality check failed for the " << configFile.getFilename() << " file."  << endl;
		}
	}
	else if(fileTypeStr==KlevanskyB0Integral3DCutoffConfigKeys::CalculationType::evaluateIntegral)
	{	
		//Check if file is written correctly
		const KlevanskyB0Integral3DCutoffFileParser config(configFile);
		if(config.validateFileQualityEvaluateIntegral())
		{	
			config.evaluateKlevanskyB0Integral3DCutoff();
		}
		else
		{
			cout << "The quality check failed for the " << configFile.getFilename() << " file."  << endl;
		}
	}
	else if (fileTypeStr==SU3NJL3DCutoffConfigKeys::CalculationType::FixedTempRhoBEqualChemPot_evaluateFirstOrderLine)
	{
		//Check if file is written correctly
		const SU3NJL3DCutoffFixedTempRhoBEqualChemPotFileParser config(configFile);
		if(config.validateFileQualityEvaluateFirstOrderLine())
		{	
			config.evaluateFirstOrderLine();
		}
		else
		{
			cout << "The quality check failed for the " << configFile.getFilename() << " file."  << endl;
		}
	}
	else
	{
		cout << "The file " << configFile.getFilename() << " does not match any known configuration! Check the FileDetails.\n";
	}
}

#include "command_line_processor.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCalculator.h"

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


	if(fileTypeStr=="evaluateSU3NJL3DCutoffVacuumMasses")
	{	
		//Check if file is written conrrectly
		const SU3NJL3DCutoffVacuumFileParser config(configFile);
		bool fileIsNice = config.validateFileQuality();
		if( fileIsNice==true )
		{	
			//evaluateSU3NJL3DCutoffVacuumMasses(configFile);
			config.evaluateVacuumMasses();
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

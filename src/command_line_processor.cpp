#include "command_line_processor.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "njl_model/n_fermion_line_integrals/KlevanskyB0Integral3DCutoffFileParser.h"

using namespace std;

int commandLineArgsProcessor(int argc, char* argv[])
{	
	string use_file_command = "use-config-file";

	// Handle command line input
    if (argc > 1) 
	{
        string command = argv[1];
        // Check for the commands
        if (command == use_file_command) 
		{
            if (argc == 3) 
			{	
				// Get name of the file from the third argument
                string configFileName = argv[2];

				// Open configuration file and parse it
				cout << "\nFeeding IniFileParser with file " << configFileName << "..." << endl;
				IniFileParser configFile(configFileName);
				selectPathBasedOnFileDetails(configFile);
            } 
			else 
			{
                cerr << "Error: No file provided after the command: " << command << endl;
                return 1;
            }
        } 
		else 
		{
            cerr << "Unknown command: " << command << endl;
            return 1;
        }
    }

	return 0;
}

void selectPathBasedOnFileDetails(const IniFileParser& configFile)
{	
    // Get the file type
    string type = configFile.getValue("FileDetails", "type");
	
	cout << "\nFileDetails:" << endl;
	cout << "type = " << type << endl;

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
	else
	{
		cout << "The file " << configFile.getFilename() << " does not match any known configuration! Check the FileDetails.\n";
	}
}

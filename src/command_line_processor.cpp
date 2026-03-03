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

void printQualityCheckFailedMessage(string filename)
{
	cout << "The quality check failed for the " << filename << " file."  << endl;
}

void selectPathBasedOnFileDetails(const IniFileParser& configFile)
{	
    // Get the file type
    string fileTypeStr = configFile.getValue("FileDetails", "type");
	
	cout << "\nFileDetails:" << endl;
	cout << "type = " << fileTypeStr << endl;

	// Check if file is written correctly and then make calculation
	if(fileTypeStr==SU3NJL3DCutoffFileParser::Vacuum::VacuumMasses::vacuumMasses)
	{	
		const SU3NJL3DCutoffFileParser::Vacuum::VacuumMasses config(configFile);
		if(config.validateFile())
		{	
			config.evaluate();
		}
		else
		{ 
			printQualityCheckFailedMessage(configFile.getFilename());
		}
	}
	else if(fileTypeStr==KlevanskyB0Integral3DCutoffFileParser::klevanskyB0Integral3DCutoff)
	{	
		const KlevanskyB0Integral3DCutoffFileParser config(configFile);
		if(config.validateFileQualityEvaluateIntegral())
		{	
			config.evaluateKlevanskyB0Integral3DCutoff();
		}
		else
		{
			printQualityCheckFailedMessage(configFile.getFilename());
		}
	}
	else if (fileTypeStr==SU3NJL3DCutoffFileParser::FixedTempRhoBEqualChemPot::FirstOrderLine::firstOrderLine)
	{
		const SU3NJL3DCutoffFileParser::FixedTempRhoBEqualChemPot::FirstOrderLine config(configFile);
		if(config.validateFile())
		{
			config.evaluate();
		}
		else
		{
			printQualityCheckFailedMessage(configFile.getFilename());
		}
	}
	else if (fileTypeStr==SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricCrossSections::isospinSymmetricCrossSections)
	{
		const SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricCrossSections config(configFile);
		if(config.validateFile())
		{	
			config.evaluate();
		}
		else
		{
			printQualityCheckFailedMessage(configFile.getFilename());
		}
	}
	else if (fileTypeStr==SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricIntegratedCrossSectionsZeroChemPot::zeroChemicalPotential )
	{	
		const SU3NJL3DCutoffFileParser::FixedChemPotTemp::IsospinSymmetricIntegratedCrossSectionsZeroChemPot config(configFile);
		if(config.validateFile())
		{	
			config.evaluate();
		}
		else
		{
			printQualityCheckFailedMessage(configFile.getFilename());
		}
	}
	else if (fileTypeStr==SU3NJL3DCutoffFileParser::FixedChemPotTemp::InMediumMassesAndThermodynamics::fileType )
	{	
		const SU3NJL3DCutoffFileParser::FixedChemPotTemp::InMediumMassesAndThermodynamics config(configFile);
		if(config.validateFile())
		{	
			config.evaluate();
		}
		else
		{
			printQualityCheckFailedMessage(configFile.getFilename());
		}
	}
	else
	{
		cout << "The file " << configFile.getFilename() << " does not match any known configuration! Check the FileDetails.\n";
	}
}

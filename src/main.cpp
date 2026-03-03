// Renan Câmara Pereira 2022-20??

#include <omp.h>
#include "command_line_processor.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedTempRhoBEqualChemPot.h"
#include "physics_utils/physical_constants.h"

#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCrossSections.h"

#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.h"

using namespace std;


int main(int argc, char* argv[])
{
	//START COUNTING TIME
    double start_s = omp_get_wtime();


    // Call the function and check its return value
	int processorResult = commandLineArgsProcessor(argc, argv);
    if ( processorResult==1 ) 
	{
        // If it returns 1, exit with failure code
        return 1;
    }
	else
	{
    	// Continue the main program execution if the function returns 0
    	std::cout << "\nCommands processed successfully, continuing execution..." << std::endl;
	}

    
    /*
    //parameter set A (Klevansky parameter set)
    double cutoff = 0.6023;
    double gs = 10.116734156126128;
    double kappa = -155.93878816540243;
    double m0u = 0.0055;
    double m0d = 0.0055;
    double m0s = 0.1407;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setA");


    double precisionVacuum = 1E-8;
    string methodVacuum = "HYBRIDS";
    double upQuarkMassGuess = 0.3;
    double strangeQuarkMassGuess = 0.5;
    double precisionVacToFinTemp = 1E-8;
    string methodVacToFinTemp = "HYBRIDS";

    double nearVacuumTemperature = 1E-5;
    double minimumTemperature = 0.120;
    double maximumTemperature = 0.300;
    int numberOfPointsFromVacToMinTemp = 200; 
    int numberOfPointsFromMinToMaxTemp = 11;
    bool largeAngleScatteringContribution = false;
    IntegratedCrossSectionApproximationMethod approximationMethod = COMPLETE_COV;
	double propagatorIntegralPrecision = 1E-5;
	double crossSectionIntegralPrecision = 1E-4;
	double integratedCrossSectionIntegralPrecision_dXdY = 1E-10;
	double integratedCrossSectionIntegralPrecision_dX = 1E-3;
    int numberOfThreads = 15;

    evaluateIsospinSymmetricIntegratedCrossSectionsWithZeroChemicalPotential(
        parameters, 
		precisionVacuum, 
		stringToMultiRootFindingMethod(methodVacuum), 
		upQuarkMassGuess, 
		strangeQuarkMassGuess, 
        precisionVacToFinTemp,
        stringToMultiRootFindingMethod(methodVacToFinTemp), 
        nearVacuumTemperature,
        minimumTemperature, 
        maximumTemperature, 
        numberOfPointsFromVacToMinTemp, 
        numberOfPointsFromMinToMaxTemp, 
        largeAngleScatteringContribution, 
        approximationMethod,
        propagatorIntegralPrecision,
        crossSectionIntegralPrecision,
        integratedCrossSectionIntegralPrecision_dXdY,
        integratedCrossSectionIntegralPrecision_dX,
        numberOfThreads
    );
*/
/*
    SU3NJL3DCutoffVacuum vacuum = SU3NJL3DCutoffVacuum::calculateVacuumMasses(
        parameters,                                    
        precisionVacuum,                                    
        methodVacuum,                                    
        mUGuess, 
        mDGuess, 
        mSGuess
    );

    double T = 0.250;
    double effChemPotU = 0.0;
    double effChemPotD = 0.0;
    double effChemPotS = 0.0;

    //find solution in the at finite temperature and fixed chemical potentials
    SU3NJL3DCutoffFixedChemPotTemp inMedium(parameters, T, effChemPotU, effChemPotD, effChemPotS);
    inMedium.solve(1E-8, HYBRIDS, vacuum.getUpQuarkEffectiveMass(), 
                                  vacuum.getDownQuarkEffectiveMass(), 
                                  vacuum.getStrangeQuarkEffectiveMass());

    cout << "inMediumSolution=" << inMedium.testSolution(1E-8) << "\n";
    cout << "Mu=" << inMedium.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << inMedium.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << inMedium.getStrangeQuarkEffectiveMass() << "GeV" << "\n";

    
    cout << "\n\n";

    double effMassU = inMedium.getUpQuarkEffectiveMass();
    double effMassD = inMedium.getDownQuarkEffectiveMass();
    double effMassS = inMedium.getStrangeQuarkEffectiveMass();


    scatteringProcess process = UUUU;
    evaluateCrossSectionProcess12To34ToFile(parameters, T, 
                                            effChemPotU, effChemPotD, effChemPotS, 
                                            effMassU, effMassD, effMassS, 
                                            1E-8, process,  
                                            false, 1E-4,
                                            20);
*/
/*
    evaluateCrossSectionsKlevanskyPaper(parameters, T, 
                                        effChemPotU, effChemPotD, effChemPotS, 
                                        effMassU, effMassD, effMassS, 
                                        1E-8,
                                        false, 1E-4,
                                        200, 14);
*/
/*
    double effMassU, effMassD, effMassS;

    double T = 0.250;
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSol = solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, T, 100, 1E-8, HYBRIDS);

    effMassU = finiteTempSol[int(finiteTempSol.size()-1)].getUpQuarkEffectiveMass();
    effMassD = finiteTempSol[int(finiteTempSol.size()-1)].getDownQuarkEffectiveMass();
    effMassS = finiteTempSol[int(finiteTempSol.size()-1)].getStrangeQuarkEffectiveMass();

    cout << "Mu=" << effMassU << "GeV" << "\t" 
         << "Md=" << effMassD << "GeV" << "\t" 
         << "Ms=" << effMassS << "GeV" << "\n";


    double chemPot = 0.100;
    vector<SU3NJL3DCutoffFixedChemPotTemp> inMediumSol = solveFromFiniteTemperatureToFiniteChemicalPotential(finiteTempSol[int(finiteTempSol.size()-1)], chemPot, 100, 1E-8, HYBRIDS);

    effMassU = inMediumSol[int(inMediumSol.size()-1)].getUpQuarkEffectiveMass();
    effMassD = inMediumSol[int(inMediumSol.size()-1)].getDownQuarkEffectiveMass();
    effMassS = inMediumSol[int(inMediumSol.size()-1)].getStrangeQuarkEffectiveMass();

    cout << "Mu=" << effMassU << "GeV" << "\t" 
         << "Md=" << effMassD << "GeV" << "\t" 
         << "Ms=" << effMassS << "GeV" << "\n";
*/


/*
    double chemPot = 0.318434158842783;
    double minimumTemperature = 0.040;
    double maximumTemperature = 0.300;
    int numberOfPointsFromVacToMinTemp = 200; 
    int numberOfPointsMinTempToChemPot = 200;
    int numberOfPointsFromMinToMaxTemp = 261;
    bool largeAngleScatteringContribution = false;
    IntegratedCrossSectionApproximationMethod approximationMethod = COMPLETE_COV;
    double propagatorIntegralPrecision = 1E-7;
    double crossSectionIntegralPrecision = 1E-4;
    double integratedCrossSectionIntegralPrecision_dXdY = 1E-12;
    double integratedCrossSectionIntegralPrecision_dX = 1E-3;
    int numberOfThreads = 14;
    evaluateIntegratedCrossSectionsWithFixedChemicalPotential(
        vacuum, 
        1E-8, 
        chemPot, 
        1E-5, 
        minimumTemperature, 
        maximumTemperature, 
        numberOfPointsFromVacToMinTemp, 
        numberOfPointsMinTempToChemPot, 
        numberOfPointsFromMinToMaxTemp, 
        largeAngleScatteringContribution, 
        approximationMethod, 
        propagatorIntegralPrecision, 
        crossSectionIntegralPrecision, 
        integratedCrossSectionIntegralPrecision_dXdY, 
        integratedCrossSectionIntegralPrecision_dX, 
        numberOfThreads
    );

    largeAngleScatteringContribution = true;
    evaluateIntegratedCrossSectionsWithFixedChemicalPotential(
        vacuum, 
        1E-8, 
        chemPot, 
        1E-5,
        minimumTemperature, 
        maximumTemperature, 
        numberOfPointsFromVacToMinTemp, 
        numberOfPointsMinTempToChemPot, 
        numberOfPointsFromMinToMaxTemp, 
        largeAngleScatteringContribution, 
        approximationMethod, 
        propagatorIntegralPrecision, 
        crossSectionIntegralPrecision, 
        integratedCrossSectionIntegralPrecision_dXdY, 
        integratedCrossSectionIntegralPrecision_dX, 
        numberOfThreads
    );
*/
/*
    double minimumTemperature = 0.120;
    double maximumTemperature = 0.300;
    int numberOfPointsFromVacToMinTemp = 200; 
    int numberOfPointsFromMinToMaxTemp = 181;
    bool largeAngleScatteringContribution = false;
    IntegratedCrossSectionApproximationMethod approximationMethod = COMPLETE_COV;
    double propagatorIntegralPrecision = 1E-7;
    double crossSectionIntegralPrecision = 1E-4;
    double integratedCrossSectionIntegralPrecision_dXdY = 1E-12;
    double integratedCrossSectionIntegralPrecision_dX = 1E-3;
    evaluateIntegratedCrossSectionsWithZeroChemicalPotential(vacuum,
                                                             1E-8,
                                                             minimumTemperature, 
                                                             maximumTemperature, 
                                                             numberOfPointsFromVacToMinTemp, 
                                                             numberOfPointsFromMinToMaxTemp, 
                                                             largeAngleScatteringContribution, 
                                                             approximationMethod,
                                                             propagatorIntegralPrecision,
                                                             crossSectionIntegralPrecision,
                                                             integratedCrossSectionIntegralPrecision_dXdY,
                                                             integratedCrossSectionIntegralPrecision_dX);

    largeAngleScatteringContribution = true;
    evaluateIntegratedCrossSectionsWithZeroChemicalPotential(vacuum,
                                                             1E-8,
                                                             minimumTemperature, 
                                                             maximumTemperature, 
                                                             numberOfPointsFromVacToMinTemp, 
                                                             numberOfPointsFromMinToMaxTemp, 
                                                             largeAngleScatteringContribution, 
                                                             approximationMethod,
                                                             propagatorIntegralPrecision,
                                                             crossSectionIntegralPrecision,
                                                             integratedCrossSectionIntegralPrecision_dXdY,
                                                             integratedCrossSectionIntegralPrecision_dX);
*/




/*
/////////////////////////////////////////////////////////////////////////////////////
//Neutron star equation of state stuff


    //parameter set (Renan Master thesis parameter set)
    double cutoff = 0.6023;
    double gs = 2*( 1.835/pow(cutoff,2) );
    double kappa = -12.360/pow(cutoff,5);
    double m0u = 0.0055;
    double m0d = 0.0055;
    double m0s = 0.1407;

    //double gOmega1 = 0.2*pow(0.5*gs, 1);
    //double gOmega2 = 5.0*pow(0.5*gs, 4);
    //double gOmega3 = -10.0*pow(0.5*gs, 7);
    //double gOmega4 = 7.0*pow(0.5*gs, 10);

    //Fix Lagrangian dimensionful couplings
    //NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_VP4Q_VP8Q_VP12Q_VP16Q, gs, kappa, gOmega1, gOmega2, gOmega3, gOmega4);

    vector<double> gOmegaAdimensional = {0.8, 1.0, -3.0, 3.0, -1.0};
    vector<double> gOmegaDimensionful = multiQuarkVPCouplingWithDimensions(gOmegaAdimensional, 0.5*gs);
    
    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_VPMULTIQ, gs, kappa, gOmegaDimensionful);


    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("renanMasterThesis");


    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, HYBRIDS, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(gapPrecision) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";


    double rhoi = 1E-5*pow(PhysicalConstants::hc_GeVfm, 3);
    double rhof = 3.80*pow(PhysicalConstants::hc_GeVfm, 3);
    int NrhoB = 8000;
    writeBetaEquilibriumEOSAtZeroTemperatureToFile(vacuum, rhoi, rhof, NrhoB, gapPrecision, HYBRIDS, "eos.dat");


    std::ofstream file;
    file.open("gOmega.dat", std::fstream::in | std::ofstream::out | std::ios::trunc);
    string gOmegaId;
    for (int i = 0; i < int( gOmegaAdimensional.size() ); ++i)
    {
        gOmegaId = gOmegaId + "_" + to_string(gOmegaAdimensional[i]);
    }
    file << gOmegaId;
    file.close();
*/


    //STOP CLOCK AND PRINT RUN TIME
    double stop_s = omp_get_wtime();
    double run_time = (stop_s-start_s);
    std::cout << "Run Time: " << run_time << std::endl;

	return 0;
}




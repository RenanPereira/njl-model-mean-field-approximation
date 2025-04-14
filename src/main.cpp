// Renan Câmara Pereira 2022-20??

#include <omp.h>
#include "command_line_processor.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedTempRhoBEqualChemPot.h"
#include "physics_utils/physical_constants.h"

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


    //////////////////////////////////////////////////////////////////////////////////////////
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

    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, DNEWTON, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";

    double minimumBaryonDensity = 1E-4*pow(hc_GeVfm,3);
    double maximumBaryonDensity = 2.00*pow(hc_GeVfm,3);
    int numberOfPoints = 2000;
    
    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint> firtOrderLine = 
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateFirstOrderLine(
        vacuum, 
        minimumBaryonDensity, 
        maximumBaryonDensity, 
        numberOfPoints, 
        1E-8, 
        DNEWTON, 
        true,
        1E-8, 
        DNEWTON, 
        0.0001, 
        1E-8
    );

    SU3NJL3DCutoffFixedTempRhoBEqualChemPot::writeFirstOrderLineToFile(vacuum, firtOrderLine, "test.dat", true);

/*
    //////////////////////////////////////////////////////////////////////////////////////////
    //parameter set B
    double cutoff = 0.586967971572559;
    double gs = 3.0257845537123/pow(cutoff, 2);
    double kappa = -11.7520217701553/pow(cutoff, 5);
    double g1 = 36.5091272818204/pow(cutoff, 8);
    double g2 = -24.3337469126998/pow(cutoff, 8);
    double m0u = 0.00600659282357854;
    double m0d = 0.00600659282357854;
    double m0s = 0.138868437136382;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_SP8Q, gs, kappa, g1, g2);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s, "setB");

    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, HYBRIDS, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";

    double minimumBaryonDensity = 1E-4*pow(hc_GeVfm,3);
    double maximumBaryonDensity = 2.00*pow(hc_GeVfm,3);
    int numberOfPoints = 2000;
    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> test =
    solveFromVacuumToFiniteBaryonDensity(vacuum, 
                                         minimumBaryonDensity, 
                                         maximumBaryonDensity, 
                                         numberOfPoints, 
                                         1E-8, HYBRIDS);
*/

/*
    //////////////////////////////////////////////////////////////////////////////////////////
    //parameter set C
    double cutoff = 0.586967971572559;
    double gs = 2.7596253718366/pow(cutoff, 2);
    double kappa = -11.7520217701553/pow(cutoff, 5);
    double g1 = 44.861215303826/pow(cutoff, 8);
    double g2 = -24.3337469126998/pow(cutoff, 8);
    double m0u = 0.00600659282357854;
    double m0d = 0.00600659282357854;
    double m0s = 0.138868437136382;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_SP8Q, gs, kappa, g1, g2);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s, "setC");

    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, HYBRIDS, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";

    double minimumBaryonDensity = 1E-4*pow(hc_GeVfm,3);
    double maximumBaryonDensity = 2.00*pow(hc_GeVfm,3);
    int numberOfPoints = 2000;
    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> test =
    solveFromVacuumToFiniteBaryonDensity(vacuum, 
                                         minimumBaryonDensity, 
                                         maximumBaryonDensity, 
                                         numberOfPoints, 
                                         1E-8, HYBRIDS);
*/

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
*/
/*
    //parameter set B
    double cutoff = 0.586967971572559;
    double gs = 3.0257845537123/pow(cutoff, 2);
    double kappa = -11.7520217701553/pow(cutoff, 5);
    double g1 = 36.5091272818204/pow(cutoff, 8);
    double g2 = -24.3337469126998/pow(cutoff, 8);
    double m0u = 0.00600659282357854;
    double m0d = 0.00600659282357854;
    double m0s = 0.138868437136382;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_SP8Q, gs, kappa, g1, g2);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setB");
*/

/*
    //parameter set C
    double cutoff = 0.586967971572559;
    double gs = 2.7596253718366/pow(cutoff, 2);
    double kappa = -11.7520217701553/pow(cutoff, 5);
    double g1 = 44.861215303826/pow(cutoff, 8);
    double g2 = -24.3337469126998/pow(cutoff, 8);
    double m0u = 0.00600659282357854;
    double m0d = 0.00600659282357854;
    double m0s = 0.138868437136382;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_SP8Q, gs, kappa, g1, g2);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setC");


    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, HYBRIDS, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";
*/

/*
    double chemPot = 0.318434158842783;
    double minimumTemperature = 0.040;
    double maximumTemperature = 0.300;
    int numberOfPointsFromVacToMinTemp = 200; 
    int numberOfPointsMinTempToChemPot = 200;
    int numberOfPointsFromMinToMaxTemp = 261;
    bool largeAngleScatteringContribution = false;
    IntegratedCrossSectionApproximationMethod approximationMethod = completeCOV;
    double propagatorIntegralPrecision = 1E-7;
    double crossSectionIntegralPrecision = 1E-4;
    double integratedCrossSectionIntegralPrecision_dXdY = 1E-12;
    double integratedCrossSectionIntegralPrecision_dX = 1E-3;
    evaluateIntegratedCrossSectionsWithFixedChemicalPotential(vacuum,
                                                              1E-8,
                                                              chemPot,
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
                                                              integratedCrossSectionIntegralPrecision_dX);

    largeAngleScatteringContribution = true;
    evaluateIntegratedCrossSectionsWithFixedChemicalPotential(vacuum,
                                                              1E-8,
                                                              chemPot,
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
                                                              integratedCrossSectionIntegralPrecision_dX);
*/
/*
    double minimumTemperature = 0.120;
    double maximumTemperature = 0.300;
    int numberOfPointsFromVacToMinTemp = 200; 
    int numberOfPointsFromMinToMaxTemp = 181;
    bool largeAngleScatteringContribution = false;
    IntegratedCrossSectionApproximationMethod approximationMethod = completeCOV;
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


    double rhoi = 1E-5*pow(hc_GeVfm, 3);
    double rhof = 3.80*pow(hc_GeVfm, 3);
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




//Renan Câmara Pereira 2022-2023

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <omp.h>
#include <algorithm>
#include "rootSolverGSL.h"
#include "Integration1DimGSL.h"
#include "generalPhysicsAndMath.h"
#include "OneFermionLineIntegral.h"
#include "TwoFermionLineIntegral.h"
#include "NJLDimensionfulCouplings.h"
#include "SU3NJL3DCutoff.h"
#include "SU3NJL3DCutoffVacuum.h"
#include "SU3NJL3DCutoffFixedChemPotTemp.h"
#include "SU3NJL3DCutoffEqualChemPotFixedTempRhoB.h"
#include "SU3NJL3DCutoffBetaEqFixedTempRhoB.h"
#include "SU3NJL3DCutoffMesonPropagators.h"
#include "SU3NJL3DCutoffDifferentialCrossSections.h"
#include "SU3NJL3DCutoffCrossSections.h"
#include "SU3NJL3DCutoffIntegratedCrossSections.h"
#include "UnitaryGroup3Dimensions.h"
#include <gsl/gsl_complex.h>


using namespace std;


int main(void)
{	
	//START COUNTING TIME
    double start_s = omp_get_wtime();


/*
    //parameter set A (Klevansky parameter set)
    double cutoff = 0.6023;
    double gs = 10.116734156126128;
    double kappa = -155.93878816540243;
    double m0u = 0.0055;
    double m0d = 0.0055;
    double m0s = 0.1407;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(interactions_4SP_det, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(cutoffEverywhere, cutoff, couplings, m0u, m0d, m0s);
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
    NJLDimensionfulCouplings couplings(interactions_4SP_det_8SP, gs, kappa, g1, g2);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(cutoffEverywhere, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setB");


    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, hybrids, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";



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
    //evaluate integrated cross sections at zero chemical potential
    evaluateIntegratedCrossSectionsWithZeroChemicalPotentialForPaper(parameters,
                                                                     0.301, 
                                                                     0.400, 
                                                                     200, 
                                                                     100, 
                                                                     false, 
                                                                     completeCOV);
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

    double gOmega1 = 0.2*pow(0.5*gs, 1);
    double gOmega2 = 8.0*pow(0.5*gs, 4);
    double gOmega3 = -10.0*pow(0.5*gs, 7);

    
    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(interactions_4SP_det_4VP_8VP_12VP, gs, kappa, gOmega1, gOmega2, gOmega3);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(cutoffEverywhere, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("renanMasterThesis");


    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, hybrids, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(gapPrecision) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";


    double rhoi = 1E-5*pow(hc_GeVfm, 3);
    double rhof = 2.50*pow(hc_GeVfm, 3);
    int NrhoB = 5000;
    writeBetaEquilibriumEOSAtZeroTemperatureToFile(vacuum, rhoi, rhof, NrhoB, gapPrecision, hybrids);
*/



    //STOP CLOCK AND PRINT RUN TIME
    double stop_s = omp_get_wtime();
    double run_time = (stop_s-start_s);
    std::cout << "Run Time: " << run_time << std::endl;

	return 0;
}




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

/*
    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, hybrids, 0.3, 0.3, 0.5);

    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";


    //solve model at zero chemical potential up to some finite temperature
    double maximumTemperature = 0.400;
    int numberOfPoints = 400;
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTSolution = 
    solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, maximumTemperature, numberOfPoints, gapPrecision, hybrids);


    //Search for meson melting points in the range of temperatures considered above
    double mesonPropertiesPrecision = 1E-7;
    double mesonMassVacuumGuess;
    double mesonWidthVacuumGuess;
    mesonState mesonID;

    mesonMassVacuumGuess = 0.2;
    mesonWidthVacuumGuess = 0.2;
    mesonID = pionPlus;
    SU3NJL3DCutoffFixedChemPotTemp meltingPointPionPlus = nondiagonalMesonMeltingPoint(vacuum, finiteTSolution, mesonID, mesonPropertiesPrecision, hybrids, mesonMassVacuumGuess, mesonWidthVacuumGuess);
    cout << meltingPointPionPlus.getTemperature() << "\n";


    ///////////////////////////////////////////////////////////////////////////
*/
    
    double T = 0.2;
    double effChemPotU = 0.5;
    double effChemPotD = 0.5;
    double effChemPotS = 0.5;
    double effMassU = 0.4;
    double effMassD = 0.5;
    double effMassS = 0.6;
    double s = 2.1*effMassS;
    scatteringProcess process = SSBarUUBar;

/*
    double test= crossSectionProcess12To34(parameters, T, 
                                           effChemPotU, effChemPotD, effChemPotS, 
                                           effMassU, effMassD, effMassS, 
                                           s, 1E-8, process, 
                                           false, 1E-4);
    cout << test << "\n";
*/

/*
    double test2 = integratedCrossSectionProcess12To34(parameters, T, 
                                                       effChemPotU, effChemPotD, effChemPotS, 
                                                       effMassU, effMassD, effMassS, 
                                                       1E-8, process, 
                                                       false, 1E-4,
                                                       1E-12, 1E-3);
*/
/*
    double test2 = integratedCrossSectionProcess12To34Klevansky(parameters, T, 
                                                       effChemPotU, effChemPotD, effChemPotS, 
                                                       effMassU, effMassD, effMassS, 
                                                       1E-8, process, 
                                                       true, 1E-4,
                                                       1E-12, 1E-3);
*/
/*
    double test2 = integratedCrossSectionProcess12To34Zhuang(parameters, T, 
                                                       effChemPotU, effChemPotD, effChemPotS, 
                                                       effMassU, effMassD, effMassS, 
                                                       1E-8, process, 
                                                       true, 1E-4,
                                                       1E-3);  
*/    


	//STOP CLOCK AND PRINT RUN TIME
    double stop_s = omp_get_wtime();
    double run_time = (stop_s-start_s);
    std::cout << "Run Time: " << run_time << std::endl;

	return 0;
}




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


    //find solution in the vacuum
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(1E-8, hybrids, 0.3, 0.3, 0.5);

    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";




    /////////////////////////////////////////////
    double maximumTemperature = 0.350;
    int numberOfPoints = 400;
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTSolution = 
    solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, maximumTemperature, numberOfPoints, 1E-8, hybrids);


    mesonState meson = pionPlus;
    double mesonPropertiesPrecision = 1E-6;
    MultiRootFindingMethod method = hybrids;
    double mesonMassVacuumGuess = 0.2;
    double mesonWidthVacuumGuess = 0.2;

    vector<SU3NJL3DCutoffMeson> mesonFiniteT = 
    mesonPropertiesFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, finiteTSolution, meson, mesonPropertiesPrecision, method, mesonMassVacuumGuess, mesonWidthVacuumGuess);


    ///////////////////////////////////////////

    for (int i = 0; i < int(mesonFiniteT.size()); ++i)
    {
        cout << mesonFiniteT[i].getTemperature() << "\t" <<
                mesonFiniteT[i].getUpQuarkEffectiveMass() << "\t" << 
                mesonFiniteT[i].getDownQuarkEffectiveMass() << "\t" << 
                mesonFiniteT[i].getStrangeQuarkEffectiveMass() << "\t" <<
                mesonFiniteT[i].getMesonMass() << "\t" <<
                mesonFiniteT[i].getMesonWidth() << "\n";
    }

    

	//STOP CLOCK AND PRINT RUN TIME
    double stop_s = omp_get_wtime();
    double run_time = (stop_s-start_s);
    std::cout << "Run Time: " << run_time << std::endl;

	return 0;
}




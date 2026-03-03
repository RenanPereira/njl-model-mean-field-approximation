#include <omp.h>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"

using namespace std;

int main()
{
	//START COUNTING TIME
    double start_s = omp_get_wtime();

    //parameter set A (Klevansky parameter set)
    double cutoff_GeV = 0.6023;
    double upQuarkCurrentMass_GeV = 0.0055;
    double downQuarkCurrentMass_GeV = 0.0055;
    double strangeQuarkCurrentMass_GeV = 0.1407;
    double fourQuarkSPCouplingCutoff2 = 3.67;
    double determinantCouplingCutoff5 = -12.36;

    double gs = fourQuarkSPCouplingCutoff2/pow(cutoff_GeV,2);
    double kappa = determinantCouplingCutoff5/pow(cutoff_GeV,5);

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(
        CUTOFF_EVERYWHERE_WITH_CTMU, 
        cutoff_GeV, 
        couplings, 
        upQuarkCurrentMass_GeV, 
        downQuarkCurrentMass_GeV, 
        strangeQuarkCurrentMass_GeV
    );
    parameters.setParameterSetName("setA");

    double precisionVacuum = 1E-8;
    string methodVacuum = "HYBRIDS";
    double upQuarkMassGuess = 0.3;
    double downQuarkMassGuess = 0.3;
    double strangeQuarkMassGuess = 0.5;

    double nearVacuumTemperature = 1E-4;
    double temperature = 0.500;
    int numberOfPoints = 2000; 
    double precisionVacToFinTemp = 1E-8;
    string methodVacToFinTemp = "HYBRIDS";

    SU3NJL3DCutoffFixedChemPotTemp::evaluateInMediumMassesAndThermodynamics(
        parameters,                                    
        precisionVacuum,                                    
        stringToMultiRootFindingMethod(methodVacuum),                                  
        upQuarkMassGuess, 
        downQuarkMassGuess,
        strangeQuarkMassGuess,
        nearVacuumTemperature,
        temperature, 
        numberOfPoints,
        precisionVacToFinTemp,
        stringToMultiRootFindingMethod(methodVacToFinTemp)
    );

    //STOP CLOCK AND PRINT RUN TIME
    double stop_s = omp_get_wtime();
    double run_time = (stop_s-start_s);
    std::cout << "Run Time: " << run_time << std::endl;

	return 0;
}

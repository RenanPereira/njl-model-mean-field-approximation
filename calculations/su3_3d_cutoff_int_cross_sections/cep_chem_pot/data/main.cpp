#include <omp.h>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.h"

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
        CUTOFF_EVERYWHERE, 
        cutoff_GeV, 
        couplings, 
        upQuarkCurrentMass_GeV, 
        downQuarkCurrentMass_GeV, 
        strangeQuarkCurrentMass_GeV
    );
    parameters.setParameterSetName("setA");

    double precisionVacuum = 1E-8;
    string methodVacuum = "DNEWTON";
    double lightQuarkMassGuess = 0.3;
    double strangeQuarkMassGuess = 0.5;

    double nearVacuumTemperature = 1E-5;
    double minimumTemperature = 0.040;
    int numberOfPointsFromVacToMinTemp = 200; 
    double precisionVacToFinTemp = 1E-8;
    string methodVacToFinTemp = "DNEWTON";

    double chemicalPotential = 0.318434158842783;
    int numberOfPointsMinTempToChemPot = 200;
    double precisionMinTempToChemPot = 1E-8;
    string methodMinTempToChemPot = "DNEWTON";
    
    double maximumTemperature = 0.300;
    int numberOfPointsFromMinToMaxTemp = 261;
	double precisionMinToMaxTemp = 1E-8;
	string methodMinToMaxTemp = "DNEWTON";
    
    bool largeAngleScatteringContribution = false;
    IntegratedCrossSectionApproximationMethod approximationMethod = COMPLETE_COV;
	double propagatorIntegralPrecision = 1E-5;
	double crossSectionIntegralPrecision = 1E-4;
	double integratedCrossSectionIntegralPrecision_dXdY = 1E-10;
	double integratedCrossSectionIntegralPrecision_dX = 1E-3;
    int numberOfThreads = 15;

    evaluateIsospinSymmetricIntegratedCrossSectionsWithFixedChemicalPotential(
        parameters, 
        precisionVacuum, 
        stringToMultiRootFindingMethod(methodVacuum), 
        lightQuarkMassGuess, 
        strangeQuarkMassGuess,
        nearVacuumTemperature,
        minimumTemperature,
        numberOfPointsFromVacToMinTemp, 
        precisionVacToFinTemp, 
        stringToMultiRootFindingMethod(methodVacToFinTemp), 
        chemicalPotential,
        numberOfPointsMinTempToChemPot,
        precisionMinTempToChemPot, 
        stringToMultiRootFindingMethod(methodMinTempToChemPot),
        maximumTemperature, 
        numberOfPointsFromMinToMaxTemp, 
        precisionMinToMaxTemp, 
        stringToMultiRootFindingMethod(methodMinToMaxTemp),
        largeAngleScatteringContribution, 
        approximationMethod,
        propagatorIntegralPrecision,
        crossSectionIntegralPrecision,
        integratedCrossSectionIntegralPrecision_dXdY,
        integratedCrossSectionIntegralPrecision_dX,
        numberOfThreads
    );

    //STOP CLOCK AND PRINT RUN TIME
    double stop_s = omp_get_wtime();
    double run_time = (stop_s-start_s);
    std::cout << "Run Time: " << run_time << std::endl;

	return 0;
}

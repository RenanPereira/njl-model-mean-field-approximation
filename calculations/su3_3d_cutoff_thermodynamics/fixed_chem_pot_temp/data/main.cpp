#include <omp.h>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"
#include <algorithm>

using namespace std;

int main()
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
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE_WITH_CTMU, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setA");

    double precisionVacuum = 1E-8;
    string methodVacuum = "HYBRIDS";
    double upQuarkMassGuess = 0.3;
    double downQuarkMassGuess = 0.3;
    double strangeQuarkMassGuess = 0.5;

    double temperature = 0.500;
    int numberOfPoints = 2000; 
    double precisionVacToFinTemp = 1E-8;
    string methodVacToFinTemp = "HYBRIDS";

	SU3NJL3DCutoffVacuum vacuum = SU3NJL3DCutoffVacuum::calculateVacuumMasses(
        parameters,                                    
        precisionVacuum,                                    
        stringToMultiRootFindingMethod(methodVacuum),                                    
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess
    );

    //solve model at zero chemical potential up to some finite temperature
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTSolution = 
    solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
        vacuum, 
        temperature, 
        numberOfPoints, 
        precisionVacToFinTemp, 
        stringToMultiRootFindingMethod(methodVacToFinTemp)
    );

    //calculate vacuum pressure
    cout << "\nCalculating the vacuum pressure...\n";
    double pressureVac = vacuum.calculatePressure();
    //cout << "vacuumPressure[GeV^4]=" << pressureVac << "\n";

    cout << "\nCalculating thermodynamics for the found solutions...\n";
    for (int i = 0; i < int(finiteTSolution.size()); ++i)
    {   
        finiteTSolution[i].setPressure(pressureVac);
        finiteTSolution[i].setEnergyDensity(-pressureVac);
        finiteTSolution[i].setEntropyDensity();
    }

    //Add solution with vacuum solution to start of the vector
    SU3NJL3DCutoffFixedChemPotTemp auxVacuum(vacuum);
    finiteTSolution.insert(finiteTSolution.begin(), auxVacuum);

    string fileName = "SU3NJL3DCutoffFixedChemPotTemp";
    fileName = fileName + "_" + finiteTSolution[0].getParametersNJL().getParameterSetName();
    fileName = fileName + "_TMin" + to_string(finiteTSolution[0].getTemperature());
    fileName = fileName + "_TMax" + to_string(finiteTSolution[finiteTSolution.size()-1].getTemperature());
    fileName = fileName + "_CP0";
    replace( fileName.begin(), fileName.end(), '.', 'p'); 
    fileName =  fileName +".dat";
    writeSolutionsToFile(finiteTSolution, fileName);

    //STOP CLOCK AND PRINT RUN TIME
    double stop_s = omp_get_wtime();
    double run_time = (stop_s-start_s);
    std::cout << "Run Time: " << run_time << std::endl;

	return 0;
}

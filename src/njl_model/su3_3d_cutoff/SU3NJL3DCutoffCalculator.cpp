#include <iostream>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCalculator.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "physics_utils/physical_constants.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedTempRhoBEqualChemPot.h"

using namespace std;


void SU3NJL3DCutoffCalculator::evaluateVacuumMasses(
    SU3NJL3DCutoffParameters& parameters,                                    
    double gapPrecision,                                    
    MultiRootFindingMethod method,                                    
    double upQuarkMassGuess, 
    double downQuarkMassGuess, 
    double strangeQuarkMassGuess
)
{   
    SU3NJL3DCutoffVacuum vacuum = SU3NJL3DCutoffVacuum::evaluateVacuumMasses(
        parameters,                                    
        gapPrecision,                                    
        method,                                    
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess
    );

    vacuum.logVacuumSolutionToFile("SU3NJL3DCutoffVacuumMasses_" + parameters.getParameterSetName() + ".dat");
}

void SU3NJL3DCutoffCalculator::evaluateFirstOrderLine(
    SU3NJL3DCutoffParameters& parameters,                                    
    double precisionVacuum,                                    
    MultiRootFindingMethod methodVacuum,                                    
    double upQuarkMassGuess, 
    double downQuarkMassGuess, 
    double strangeQuarkMassGuess,
    double minimumBaryonDensity_fmMinus3, 
    double maximumBaryonDensity_fmMinus3, 
    int numberOfPoints,
    double precisionZeroTempSol,
    MultiRootFindingMethod methodZeroTempSol,
    bool storeZeroTempSolToFile,
    double precisionTransitionPointSol, 
    MultiRootFindingMethod methodTransitionPointSol,
    double deltaT,
    double massDifferenceCEP
)
{
    SU3NJL3DCutoffVacuum vacuum = SU3NJL3DCutoffVacuum::evaluateVacuumMasses(
        parameters,                                    
        precisionVacuum,                                    
        methodVacuum,                                    
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess
    );

    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint> firtOrderLine = 
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateFirstOrderLine(
        vacuum, 
        minimumBaryonDensity_fmMinus3*pow(hc_GeVfm,3), 
        maximumBaryonDensity_fmMinus3*pow(hc_GeVfm,3), 
        numberOfPoints, 
        precisionZeroTempSol, 
        methodZeroTempSol, 
        storeZeroTempSolToFile,
        precisionTransitionPointSol, 
        methodTransitionPointSol, 
        deltaT, 
        massDifferenceCEP
    );

    SU3NJL3DCutoffFixedTempRhoBEqualChemPot::writeToFile(
        vacuum, 
        firtOrderLine, 
        "SU3NJL3DCutoffFirstOrderLine_" + parameters.getParameterSetName() + ".dat",
        true
    );
}

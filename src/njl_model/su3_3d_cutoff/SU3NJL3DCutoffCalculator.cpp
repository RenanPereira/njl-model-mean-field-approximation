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

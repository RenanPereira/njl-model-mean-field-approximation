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


void SU3NJL3DCutoffCalculator::evaluateVacuumMasses(const IniFileParser& config)
{   
    // Get SU3 NJL 3D Cutoff Model Parameters
    cout << "\nSU3NJL3DCutoffModelParameters:" << endl;

    string parameterSetName = config.getValue("SU3NJL3DCutoffModelParameters", "parameterSetName");
    string regularizationScheme = config.getValue("SU3NJL3DCutoffModelParameters", "regularizationScheme");
    double cutoffInGeV = config.getDouble("SU3NJL3DCutoffModelParameters", "cutoffInGeV");
    double upQuarkCurrentMassInGeV = config.getDouble("SU3NJL3DCutoffModelParameters", "upQuarkCurrentMassInGeV");
    double downQuarkCurrentMassInGeV = config.getDouble("SU3NJL3DCutoffModelParameters", "downQuarkCurrentMassInGeV");
    double strangeQuarkCurrentMassInGeV = config.getDouble("SU3NJL3DCutoffModelParameters", "strangeQuarkCurrentMassInGeV");
    
    cout << "parameterSetName = " << parameterSetName << endl;
    cout << "regularizationScheme = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << "cutoffInGeV = " << cutoffInGeV << endl;
    cout << "upQuarkCurrentMassInGeV = " << upQuarkCurrentMassInGeV << endl;
    cout << "downQuarkCurrentMassInGeV = " << downQuarkCurrentMassInGeV << endl;
    cout << "strangeQuarkCurrentMassInGeV = " << strangeQuarkCurrentMassInGeV << endl;

    // Get SU3 NJL 3D Cutoff Dimensionful Couplings
    cout << "\nSU3NJL3DCutoffGapEquationsVacuumParameters:" << endl;
    NJLDimensionfulCouplings couplings = SU3NJL3DCutoffFileParser::extractDimensionfulCouplings(config);

    // SU3 NJL 3D Cutoff Gap Equations Vacuum Parameters
    cout << "\nSU3NJL3DCutoffGapEquationsVacuumParameters: " << endl;

    double gapPrecision = config.getDouble("SU3NJL3DCutoffGapEquationsVacuumParameters", "gapPrecision");
    string rootFindingMethod = config.getValue("SU3NJL3DCutoffGapEquationsVacuumParameters", "rootFindingMethod");
    double upQuarkMassGuess = config.getDouble("SU3NJL3DCutoffGapEquationsVacuumParameters", "upQuarkMassGuess");
    double downQuarkMassGuess = config.getDouble("SU3NJL3DCutoffGapEquationsVacuumParameters", "downQuarkMassGuess");
    double strangeQuarkMassGuess = config.getDouble("SU3NJL3DCutoffGapEquationsVacuumParameters", "strangeQuarkMassGuess");

    cout << "gapPrecision = " << gapPrecision << endl;
    cout << "rootFindingMethod = " << toStringMultiRootFindingMethod(stringToMultiRootFindingMethod(rootFindingMethod)) << endl;
    cout << "upQuarkMassGuess = " << upQuarkMassGuess << endl;
    cout << "downQuarkMassGuess = " << downQuarkMassGuess << endl;
    cout << "strangeQuarkMassGuess = " << strangeQuarkMassGuess << endl;

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(
        stringToNJL3DCutoffRegularizationScheme(regularizationScheme),                                 
        cutoffInGeV,                                 
        couplings,                                 
        upQuarkCurrentMassInGeV,                                 
        downQuarkCurrentMassInGeV,                                 
        strangeQuarkCurrentMassInGeV
    );
    parameters.setParameterSetName(parameterSetName);

    // Solve model in the vacuum
    evaluateVacuumMasses(
        parameters,                                
        gapPrecision,                                
        stringToMultiRootFindingMethod(rootFindingMethod),                                
        upQuarkMassGuess,                                
        downQuarkMassGuess,                                
        strangeQuarkMassGuess
    );
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

#include <iostream>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCalculator.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"

using namespace std;


void evaluateSU3NJL3DCutoffVacuumMasses(SU3NJL3DCutoffParameters& parameters,
                                        double gapPrecision,
                                        MultiRootFindingMethod method,
                                        double upQuarkMassGuess, double downQuarkMassGuess, double strangeQuarkMassGuess)
{   
    // Solve model in the vacuum
    cout << "\nSolving the SU3 NJL model, regularized by a 3D Cutoff, in vacuum...\n";

    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, 
                 method, 
                 upQuarkMassGuess, 
                 downQuarkMassGuess, 
                 strangeQuarkMassGuess);

    double Mu = vacuum.getUpQuarkEffectiveMass();
    double Md = vacuum.getDownQuarkEffectiveMass();
    double Ms = vacuum.getStrangeQuarkEffectiveMass();

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(gapPrecision) << "\n";
    cout << "Mu[GeV] = " << Mu << "\n" 
         << "Md[GeV] = " << Md << "\n" 
         << "Ms[GeV] = " << Ms << "\n";

    vacuum.logVacuumSolutionToFile("SU3NJL3DCutoffVacuumMasses_" + parameters.getParameterSetName() + ".dat");
}


void evaluateSU3NJL3DCutoffVacuumMasses(const IniFileParser& config)
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
    NJLDimensionfulCouplings couplings = SU3NJL3DCutoffVacuumFileParser::extractDimensionfulCouplings(config);


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
    SU3NJL3DCutoffParameters parameters(stringToNJL3DCutoffRegularizationScheme(regularizationScheme), 
                                        cutoffInGeV, 
                                        couplings, 
                                        upQuarkCurrentMassInGeV, 
                                        downQuarkCurrentMassInGeV, 
                                        strangeQuarkCurrentMassInGeV);
    parameters.setParameterSetName(parameterSetName);


    // Solve model in the vacuum
    evaluateSU3NJL3DCutoffVacuumMasses(parameters, 
                                       gapPrecision, 
                                       stringToMultiRootFindingMethod(rootFindingMethod), 
                                       upQuarkMassGuess, 
                                       downQuarkMassGuess, 
                                       strangeQuarkMassGuess);
}

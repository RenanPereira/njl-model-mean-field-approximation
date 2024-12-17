#include <iostream>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCalculator.h"
#include "SU3NJL3DCutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"

using namespace std;


void evaluateSU3NJL3DCutoffVacuumMasses(const IniFileParser& config)
{   
    // Get SU3 NJL 3D Cutoff Model Parameters
    cout << "\nSU3NJL3DCutoffModelParameters:" << endl;

    string parameterSetName = config.getValue("SU3NJL3DCutoffModelParameters", "parameterSetName");
    string regularizationScheme = config.getValue("SU3NJL3DCutoffModelParameters", "NJL3DCutoffRegularizationScheme");
    double cutoff_GeV = config.getDouble("SU3NJL3DCutoffModelParameters", "cutoff_GeV");
    double upCurrentQuarkMass_GeV = config.getDouble("SU3NJL3DCutoffModelParameters", "upCurrentQuarkMass_GeV");
    double downCurrentQuarkMass_GeV = config.getDouble("SU3NJL3DCutoffModelParameters", "downCurrentQuarkMass_GeV");
    double strangeCurrentQuarkMass_GeV = config.getDouble("SU3NJL3DCutoffModelParameters", "strangeCurrentQuarkMass_GeV");
    
    cout << "parameterSetName = " << parameterSetName << endl;
    cout << "NJL3DCutoffRegularizationScheme = " << toStringNJL3DCutoffRegularizationScheme(stringToNJL3DCutoffRegularizationScheme(regularizationScheme)) << endl;
    cout << "cutoff_GeV = " << cutoff_GeV << endl;
    cout << "upCurrentQuarkMass_GeV = " << upCurrentQuarkMass_GeV << endl;
    cout << "downCurrentQuarkMass_GeV = " << downCurrentQuarkMass_GeV << endl;
    cout << "strangeCurrentQuarkMass_GeV = " << strangeCurrentQuarkMass_GeV << endl;


    // Get SU3 NJL 3D Cutoff Dimensionful Couplings
    NJLDimensionfulCouplings couplings = extractSU3NJL3DCutoffDimensionfulCouplings(config);


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
                                        cutoff_GeV, 
                                        couplings, 
                                        upCurrentQuarkMass_GeV, 
                                        downCurrentQuarkMass_GeV, 
                                        strangeCurrentQuarkMass_GeV);
    parameters.setParameterSetName(parameterSetName);


    // Solve model in the vacuum
    cout << "\nSolving the SU3 NJL model, regularized by a 3D Cutoff, in vacuum...\n";

    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, 
                 stringToMultiRootFindingMethod(rootFindingMethod), 
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

    vacuum.logVacuumSolutionToFile("SU3NJL3DCutoffVacuumMasses_" + parameterSetName+".dat");
}
#include <iostream>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCalculator.h"
#include "SU3NJL3DCutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"

using namespace std;


void evaluateSU3NJL3DCutoffVacuumMasses(const IniFileParser& config)
{   
    // SU3NJL3DCutoffModelParameters
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


    //////////////////////////////////////////

    // NJLDimensionfulCouplings
    cout << "\nNJLDimensionfulCouplings: " << endl;

    string interactionTerms = config.getValue("NJLDimensionfulCouplings", "lagrangianInteractions");
    cout << "interactionTerms = " << toStringLagrangianInteractions(fromStringLagrangianInteractions(interactionTerms)) << endl;

    double fourQuarkSPCouplingCutoff2 = config.getDouble("NJLDimensionfulCouplings", "fourQuarkSPCouplingCutoff2");
    cout << "fourQuarkSPCouplingCutoff2 = " << fourQuarkSPCouplingCutoff2 << endl;

    double determinantCouplingCutoff5 = config.getDouble("NJLDimensionfulCouplings", "determinantCouplingCutoff5");
    cout << "determinantCouplingCutoff5 = " << determinantCouplingCutoff5 << endl;

    double eightQuarkSPOziViolatingCouplingCutoff8 = config.getDouble("NJLDimensionfulCouplings", "eightQuarkSPOziViolatingCouplingCutoff8");
    cout << "eightQuarkSPOziViolatingCouplingCutoff8 = " << eightQuarkSPOziViolatingCouplingCutoff8 << endl;

    double eightQuarkSPNonOziViolatingCouplingCutoff8 = config.getDouble("NJLDimensionfulCouplings", "eightQuarkSPNonOziViolatingCouplingCutoff8");
    cout << "eightQuarkSPNonOziViolatingCouplingCutoff8 = " << eightQuarkSPNonOziViolatingCouplingCutoff8 << endl;

    //Fix Lagrangian couplings
    double gs = fourQuarkSPCouplingCutoff2/pow(cutoff_GeV, 2);
    double kappa = determinantCouplingCutoff5/pow(cutoff_GeV, 5);
    double g1 = eightQuarkSPOziViolatingCouplingCutoff8/pow(cutoff_GeV, 8);
    double g2 = eightQuarkSPNonOziViolatingCouplingCutoff8/pow(cutoff_GeV, 8);

    NJLDimensionfulCouplings couplings(fromStringLagrangianInteractions(interactionTerms), gs, kappa, g1, g2);
    
    //////////////////////////////////////////


    // NJLDimensionfulCouplings
    cout << "\nNJLDimensionfulCouplings: " << endl;

    double gapPrecision = config.getDouble("SU3NJL3DCutoffGapEquationsVacuum", "gapPrecision");
    string rootFindingMethod = config.getValue("SU3NJL3DCutoffGapEquationsVacuum", "rootFindingMethod");
    double upQuarkMassGuess = config.getDouble("SU3NJL3DCutoffGapEquationsVacuum", "upQuarkMassGuess");
    double downQuarkMassGuess = config.getDouble("SU3NJL3DCutoffGapEquationsVacuum", "downQuarkMassGuess");
    double strangeQuarkMassGuess = config.getDouble("SU3NJL3DCutoffGapEquationsVacuum", "strangeQuarkMassGuess");

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


    //solve model in the vacuum
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, 
                 stringToMultiRootFindingMethod(rootFindingMethod), 
                 upQuarkMassGuess, 
                 downQuarkMassGuess, 
                 strangeQuarkMassGuess);

    double Mu = vacuum.getUpQuarkEffectiveMass();
    double Md = vacuum.getDownQuarkEffectiveMass();
    double Ms = vacuum.getStrangeQuarkEffectiveMass();

    cout << "\nVacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(gapPrecision) << "\n";
    cout << "Mu[GeV] = " << Mu << "\n" 
         << "Md[GeV] = " << Md << "\n" 
         << "Ms[GeV] = " << Ms << "\n";

    vacuum.logVacuumSolutionToFile("SU3NJL3DCutoffVacuumMasses_" + parameterSetName+".dat");
}
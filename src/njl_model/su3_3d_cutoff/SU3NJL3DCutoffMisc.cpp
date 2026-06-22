#include <iostream>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffMisc.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCrossSections.h"

using std::cout;
using std::vector;


void someVacuumAndThermalPropertiesKlevanskyParameterSet()
{
    //parameter set A (Klevansky parameter set)
    double cutoff = 0.6023;
    double gs = 10.116734156126128;
    double kappa = -155.93878816540243;
    double m0u = 0.0055;
    double m0d = 0.0055;
    double m0s = 0.1407;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(LagrangianInteractions::SP4Q_DET2NFQ, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(NJL3DCutoffRegularizationScheme::CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setA");

    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, MultiRootFindingMethod::HYBRIDS, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";

    cout << "Pseudoscalar meson masses:\n";

    SU3NJL3DCutoffMeson pionPlusMassSolution = vacuum.calculateMesonMassAndWidth(MesonState::pionPlus, 1E-7, MultiRootFindingMethod::HYBRIDS, 0.2, 0.2);
    cout << pionPlusMassSolution.getMesonMass() << "\t" << pionPlusMassSolution.getMesonWidth() << "\n";

    SU3NJL3DCutoffMeson kaonPlusMassSolution = vacuum.calculateMesonMassAndWidth(MesonState::kaonPlus, 1E-7, MultiRootFindingMethod::HYBRIDS, 0.4, 0.2);
    cout << kaonPlusMassSolution.getMesonMass() << "\t" << kaonPlusMassSolution.getMesonWidth() << "\n";

    SU3NJL3DCutoffMeson diagonalMassSolution;
    diagonalMassSolution = vacuum.calculateMesonMassAndWidth(MesonState::diagonalPseudoscalars, 1E-7, MultiRootFindingMethod::HYBRIDS, 0.5, 0.2);
    cout << diagonalMassSolution.getMesonMass() << "\t" << diagonalMassSolution.getMesonWidth() << "\n";

    diagonalMassSolution = vacuum.calculateMesonMassAndWidth(MesonState::diagonalPseudoscalars, 1E-7, MultiRootFindingMethod::HYBRIDS, 1.0, 0.2);
    cout << diagonalMassSolution.getMesonMass() << "\t" << diagonalMassSolution.getMesonWidth() << "\n";

     cout << "Scalar meson masses:\n";

    SU3NJL3DCutoffMeson sigmaPionPlusMassSolution = vacuum.calculateMesonMassAndWidth(MesonState::sigmaPionPlus, 1E-7, MultiRootFindingMethod::HYBRIDS, 0.2, 0.2);
    cout << sigmaPionPlusMassSolution.getMesonMass() << "\t" << sigmaPionPlusMassSolution.getMesonWidth() << "\n";

    //solve model at zero chemical potential up to some finite temperature
    double nearVacuumTemperature = 1E-4;
    double maximumTemperature = 0.400;
    int numberOfPoints = 400;
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTSolution = 
    SU3NJL3DCutoffFixedChemPotTemp::solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
        vacuum, 
        nearVacuumTemperature,
        maximumTemperature, 
        numberOfPoints, 
        gapPrecision, 
        MultiRootFindingMethod::HYBRIDS
    );

    //Search for meson melting points in the range of temperatures considered above
    double mesonPropertiesPrecision = 1E-7;
    double mesonMassVacuumGuess;
    double mesonWidthVacuumGuess;
    MesonState mesonID;

    mesonMassVacuumGuess = 0.2;
    mesonWidthVacuumGuess = 0.2;
    mesonID = MesonState::pionPlus;
    SU3NJL3DCutoffFixedChemPotTemp meltingPointPionPlus = nondiagonalMesonMeltingPoint(
        vacuum, 
        finiteTSolution, 
        mesonID, 
        mesonPropertiesPrecision, 
        MultiRootFindingMethod::HYBRIDS, 
        mesonMassVacuumGuess, 
        mesonWidthVacuumGuess
    );
    cout << "pionPlus TMott: " << meltingPointPionPlus.getTemperature() << "\n";

    mesonMassVacuumGuess = 0.5;
    mesonWidthVacuumGuess = 0.2;
    mesonID = MesonState::kaonPlus;
    SU3NJL3DCutoffFixedChemPotTemp meltingPointKaonPlus = nondiagonalMesonMeltingPoint(
        vacuum, 
        finiteTSolution, 
        mesonID, 
        mesonPropertiesPrecision, 
        MultiRootFindingMethod::HYBRIDS, 
        mesonMassVacuumGuess, 
        mesonWidthVacuumGuess
    );
    cout << "kaonPlus TMott: " << meltingPointKaonPlus.getTemperature() << "\n";
}

//Evaluate Cross section for paper at finite chemical potential using Klevansky parameter set
void evaluateCrossSectionsPaperWithKlevanskyParameterSet(
    double nearVacuumTemperature, 
    double T, 
    double chemPot, 
    int numberOfCrossSectionPoints, 
    int numberOfThreads
)
{
    //define Klevansky parameters
    //parameter set A (Klevansky parameter set)
    double cutoff = 0.6023;
    double gs = 10.116734156126128;
    double kappa = -155.93878816540243;
    double m0u = 0.0055;
    double m0d = 0.0055;
    double m0s = 0.1407;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(LagrangianInteractions::SP4Q_DET2NFQ, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(NJL3DCutoffRegularizationScheme::CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);

    //numerical precisions
    double gapPrecision = 1E-8;
    double mesonPropagatorIntegralPrecision = 1E-8;
    double crossSectionIntegralPrecision = 1E-4;

    //find solution in the vacuum for the Klevansky parameter set
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, MultiRootFindingMethod::HYBRIDS, 0.3, 0.3, 0.5);

    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";

    //solve gap equation from the vacuum up to finite temperature
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSol = SU3NJL3DCutoffFixedChemPotTemp::solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
        vacuum, 
        nearVacuumTemperature,
        T, 
        100, 
        gapPrecision, 
        MultiRootFindingMethod::HYBRIDS
    );

    double effMassU, effMassD, effMassS;
    effMassU = finiteTempSol[int(finiteTempSol.size()-1)].getUpQuarkEffectiveMass();
    effMassD = finiteTempSol[int(finiteTempSol.size()-1)].getDownQuarkEffectiveMass();
    effMassS = finiteTempSol[int(finiteTempSol.size()-1)].getStrangeQuarkEffectiveMass();

    cout << "Mu=" << effMassU << "GeV" << "\t" 
         << "Md=" << effMassD << "GeV" << "\t" 
         << "Ms=" << effMassS << "GeV" << "\n";

    //solve gap equation from the finite temperature up to finite chemical potential
    if( chemPot>0.0 )
    {
        vector<SU3NJL3DCutoffFixedChemPotTemp> inMediumSol = SU3NJL3DCutoffFixedChemPotTemp::solveFromFiniteTemperatureToFiniteChemicalPotential(
            finiteTempSol[int(finiteTempSol.size()-1)], 
            chemPot, 
            100, 
            gapPrecision, 
            MultiRootFindingMethod::HYBRIDS
        );

        effMassU = inMediumSol[int(inMediumSol.size()-1)].getUpQuarkEffectiveMass();
        effMassD = inMediumSol[int(inMediumSol.size()-1)].getDownQuarkEffectiveMass();
        effMassS = inMediumSol[int(inMediumSol.size()-1)].getStrangeQuarkEffectiveMass();

        cout << "Mu=" << effMassU << "GeV" << "\t" 
             << "Md=" << effMassD << "GeV" << "\t" 
             << "Ms=" << effMassS << "GeV" << "\n";
    }

    //evaluate cross sections
    evaluateCrossSectionsEqualLightMassesEqualChemicalPotential(
        parameters, T, 
        chemPot, chemPot, chemPot, 
        effMassU, effMassD, effMassS, 
        mesonPropagatorIntegralPrecision,
        false, crossSectionIntegralPrecision,
        numberOfCrossSectionPoints,
        numberOfThreads
    );
}

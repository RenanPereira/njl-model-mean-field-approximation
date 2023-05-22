
/*
    double pressureVac = vacuum.calculatePressure();
    double energyVac = vacuum.calculateEnergyDensity();
    cout << "pressure=" << pressureVac << "\n";
    cout << "energyDensity=" << energyVac << "\n";

    double pressureVacElec = vacuum.calculateVacuumPressureElectrons(electronMass_GeV);
    cout << "pressureElec=" << pressureVacElec << "\n";
*/

/*
    SU3NJL3DCutoffEqualChemPotFixedTempRhoB inMedium(parameters);
    inMedium.setTemperature(0.0);
    inMedium.setBaryonDensity(0.0);

    //guesses for i=0
    double MuGuess = vacuum.getUpQuarkEffectiveMass();
    double MdGuess = vacuum.getDownQuarkEffectiveMass();
    double MsGuess = vacuum.getStrangeQuarkEffectiveMass();
    double effectiveCPGuess = vacuum.getUpQuarkEffectiveMass()*(1+1E-4);

    double rhoi = 1E-4*pow(hc_GeVfm,3);
    double rhof = 2.00*pow(hc_GeVfm,3);
    int NrhoB = 1000;
    double drhoB = (rhof-rhoi)/(NrhoB-1);
    for (int i = 0; i < NrhoB; ++i)
    {   
        //step in density
        double rho_B = rhoi + i*drhoB;
        inMedium.setBaryonDensity(rho_B);

        //find quark masses and effective chemical potential
        inMedium.solve(1E-8, hybrids, MuGuess, MdGuess, MsGuess, effectiveCPGuess);

        //guesses for next step
        MuGuess = inMedium.getUpQuarkEffectiveMass();
        MdGuess = inMedium.getDownQuarkEffectiveMass();
        MsGuess = inMedium.getStrangeQuarkEffectiveMass();
        effectiveCPGuess = inMedium.getQuarkEffectiveChemicalPotential();  

        //cout << "testSolution=" << inMedium.testSolution() << "\n";
        cout << rho_B/pow(hc_GeVfm,3) << "\t" 
             << inMedium.getUpQuarkEffectiveMass() << "\t"
             << inMedium.getDownQuarkEffectiveMass() << "\t"
             << inMedium.getStrangeQuarkEffectiveMass() << "\n";  


        double pressureMed = inMedium.calculatePressure(pressureVac);
        double energyMed = inMedium.calculateEnergyDensity(energyVac);
        double entropyMed = inMedium.calculateEntropyDensity();
        cout << "pressure=" << pressureMed << "\n";
        cout << "energyDensity=" << energyMed << "\n";
        cout << "entropyDensity=" << entropyMed << "\n";
    }
*/


/*
    vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> betaEqSolutions;
    addVacuumSolution(vacuum, electronMass_GeV, pressureVac, pressureVacElec, betaEqSolutions);


    //guesses for i=0
    double mUGuess = vacuum.getUpQuarkEffectiveMass();
    double mDGuess = vacuum.getDownQuarkEffectiveMass();
    double mSGuess = vacuum.getStrangeQuarkEffectiveMass();
    double effCPUGuess = mUGuess + 1E-3;
    double effCPDGuess = mUGuess + 1E-3;
    double effCPSGuess = mUGuess + 1E-3;

    double rhoi = 1E-5*pow(hc_GeVfm,3);
    double rhof = 2.00*pow(hc_GeVfm,3);
    int NrhoB = 5000;
    double temperature = 0.0;

    double drhoB = (rhof-rhoi)/(NrhoB-1);
    for (int i = 0; i < NrhoB; ++i)
    {   
        //step in density
        double rhoB = rhoi + i*drhoB;

        //create object
        SU3NJL3DCutoffBetaEqFixedTempRhoB betaEq(parameters, electronMass_GeV, temperature, rhoB);

        //find quark masses and effective chemical potential
        betaEq.solve(1E-8, hybrids, mUGuess, mDGuess, mSGuess, effCPUGuess, effCPDGuess, effCPSGuess);

        //guesses for next step
        mUGuess = betaEq.getUpQuarkEffectiveMass();
        mDGuess = betaEq.getDownQuarkEffectiveMass();
        mSGuess = betaEq.getStrangeQuarkEffectiveMass();
        effCPUGuess = betaEq.getUpQuarkEffectiveChemicalPotential();
        effCPDGuess = betaEq.getDownQuarkEffectiveChemicalPotential();
        effCPSGuess = betaEq.getStrangeQuarkEffectiveChemicalPotential();

        //calculate thermodynamics
        betaEq.setBetaEqThermodynamics(pressureVac, pressureVacElec);


        //print to console
        cout << rhoB/pow(hc_GeVfm,3) << "\t"
             << betaEq.getUpQuarkEffectiveMass() << "\t"
             << betaEq.getDownQuarkEffectiveMass() << "\t"
             << betaEq.getStrangeQuarkEffectiveMass() << "\t"
             << betaEq.getUpQuarkEffectiveChemicalPotential() << "\t"
             << betaEq.getDownQuarkEffectiveChemicalPotential() << "\t"
             << betaEq.getStrangeQuarkEffectiveChemicalPotential() << "\n";


        //push to solutions vector
        betaEqSolutions.push_back(betaEq);
    }


    //save all information in file
    writeSolutionsToFile(betaEqSolutions, "solutions.dat", true);


    //find chiral transition
    vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> transitionPoints = findChiralTransitionPointsFixedTemperature(betaEqSolutions, 1E-8, dnewton);


    //save EOS to file: if it has first order phase transition, save only after restoration 
    if ( int(transitionPoints.size())>0 )
    {   
        double minRhoB = transitionPoints[1].getBaryonDensity();
        writeEOSToFile(betaEqSolutions, "eos.dat", true, minRhoB);
    }
    else
    {
        writeEOSToFile(betaEqSolutions, "eos.dat", true);
    }
*/


/*
    double test2OG = integratedCrossSectionProcess12To34Zhuang(parameters, 0.7, 
                                                    0.25, 0.35, 0.45, 
                                                    0.3, 0.4, 0.5, 
                                                    1E-8, SBarSBarSBarSBar, 
                                                    false, 1E-4,
                                                    1E-5);

    cout << test2OG << "\n";
*/
/*
    int a = 6;
    int b = 1;
    int c = 5;
    gsl_complex fabc = unitaryGroup3DCalculateAntisymmetricStructureConstant(a, b, c);

    cout << GSL_REAL(fabc) << "\t" << GSL_IMAG(fabc) << "\n\n";    
*/  



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Evaluate Cross section for paper at finite chemical potential using Klevansky parameter set

void evaluateCrossSectionsPaperFiniteChemicalPotentialWithKlevanskyParameterSet()
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
    NJLDimensionfulCouplings couplings(interactions_4SP_det, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(cutoffEverywhere, cutoff, couplings, m0u, m0d, m0s);


    //find solution in the vacuum for the Klevansky parameter set
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(1E-8, hybrids, 0.3, 0.3, 0.5);

    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";

    //solve gap equation from the vacuum up to finite temperature
    double T = 0.250;
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSol = solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, T, 100, 1E-8, hybrids);

    double effMassU, effMassD, effMassS;
    effMassU = finiteTempSol[int(finiteTempSol.size()-1)].getUpQuarkEffectiveMass();
    effMassD = finiteTempSol[int(finiteTempSol.size()-1)].getDownQuarkEffectiveMass();
    effMassS = finiteTempSol[int(finiteTempSol.size()-1)].getStrangeQuarkEffectiveMass();

    cout << "Mu=" << effMassU << "GeV" << "\t" 
         << "Md=" << effMassD << "GeV" << "\t" 
         << "Ms=" << effMassS << "GeV" << "\n";


    //solve gap equation from the finite temperature up to finite chemical potential
    double chemPot = 0.000;
    vector<SU3NJL3DCutoffFixedChemPotTemp> inMediumSol = solveFromFiniteTemperatureToFiniteChemicalPotential(finiteTempSol[int(finiteTempSol.size()-1)], chemPot, 100, 1E-8, hybrids);

    effMassU = inMediumSol[int(inMediumSol.size()-1)].getUpQuarkEffectiveMass();
    effMassD = inMediumSol[int(inMediumSol.size()-1)].getDownQuarkEffectiveMass();
    effMassS = inMediumSol[int(inMediumSol.size()-1)].getStrangeQuarkEffectiveMass();

    cout << "Mu=" << effMassU << "GeV" << "\t" 
         << "Md=" << effMassD << "GeV" << "\t" 
         << "Ms=" << effMassS << "GeV" << "\n";


    //evaluate cross sections
    evaluateCrossSectionsPaperFiniteChemicalPotential(parameters, T, 
                                                      chemPot, chemPot, chemPot, 
                                                      effMassU, effMassD, effMassS, 
                                                      1E-8,
                                                      false, 1E-4,
                                                      200);
}







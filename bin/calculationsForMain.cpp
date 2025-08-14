
/*
    double pressureVac = vacuum.calculatePressure();
    double energyVac = vacuum.calculateEnergyDensity();
    cout << "pressure=" << pressureVac << "\n";
    cout << "energyDensity=" << energyVac << "\n";

    double pressureVacElec = vacuum.calculateVacuumPressureElectrons(electronMass_GeV);
    cout << "pressureElec=" << pressureVacElec << "\n";
*/

/*
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot inMedium(parameters);
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
        inMedium.solve(1E-8, HYBRIDS, MuGuess, MdGuess, MsGuess, effectiveCPGuess);

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
        betaEq.solve(1E-8, HYBRIDS, mUGuess, mDGuess, mSGuess, effCPUGuess, effCPDGuess, effCPSGuess);

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
    vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> transitionPoints = findChiralTransitionPointsFixedTemperature(betaEqSolutions, 1E-8, DNEWTON);


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



/*

//Code used to test EoS in beta equilibrium versus old code mark_6: everything checks out!

/////////////////////////////////////////////////////////////////////////////////////
//Neutron star equation of state stuff


    //parameter set (Renan Master thesis parameter set)
    double cutoff = 0.6023;
    double gs = 2*( 1.835/pow(cutoff,2) );
    double kappa = -12.360/pow(cutoff,5);
    double m0u = 0.0055;
    double m0d = 0.0055;
    double m0s = 0.1407;

    double gOmega1 = 1.0*pow(0.5*gs, 1);
    double gRho1 = 1.0*pow(0.5*gs, 1);
    double gOmega2 = 10.0*pow(0.5*gs, 4);
    double gRho2 = 10.0*pow(0.5*gs, 4);
    double gOmegaRho = 10.0*pow(0.5*gs, 4);
    double gSigmaOmega = 10.0*pow(0.5*gs, 4);
    double gSigmaRho = 5.0*pow(0.5*gs, 4);
    
    //double gOmega3 = -10.0*pow(0.5*gs, 7);

    
    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_VP4Q_VIPI4Q_VP8Q_VIPI8Q_VPVIPI8Q_SPVP8Q_SPVIPI8Q, 
                                       gs, kappa, gOmega1, gRho1, gOmega2, gRho2, gOmegaRho, gSigmaOmega, gSigmaRho);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("renanMasterThesis");


    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, HYBRIDS, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(gapPrecision) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";


    double rhoi = 1E-5*pow(hc_GeVfm, 3);
    double rhof = 2.50*pow(hc_GeVfm, 3);
    int NrhoB = 5000;
    writeBetaEquilibriumEOSAtZeroTemperatureToFile(vacuum, rhoi, rhof, NrhoB, gapPrecision, HYBRIDS);




*/




/*
// Cross sections for paper

    //parameter set A (Klevansky parameter set)
    double cutoff = 0.6023;
    double gs = 10.116734156126128;
    double kappa = -155.93878816540243;
    double m0u = 0.0055;
    double m0d = 0.0055;
    double m0s = 0.1407;

    double precisionVacuum = 1E-8;
    MultiRootFindingMethod methodVacuum = DNEWTON;
    double mUGuess = 0.3;
    double mDGuess = 0.3;
    double mSGuess = 0.5;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setA");
    
    SU3NJL3DCutoffVacuum vacuum = SU3NJL3DCutoffVacuum::calculateVacuumMasses(
        parameters,                                    
        precisionVacuum,                                    
        methodVacuum,                                    
        mUGuess, 
        mDGuess, 
        mSGuess
    );

    double T = 0.250;
    double effChemPotU = 0.0;
    double effChemPotD = 0.0;
    double effChemPotS = 0.0;

    //find solution in the at finite temperature and fixed chemical potentials
    SU3NJL3DCutoffFixedChemPotTemp inMedium(parameters, T, effChemPotU, effChemPotD, effChemPotS);
    inMedium.solve(1E-8, HYBRIDS, vacuum.getUpQuarkEffectiveMass(), 
                                  vacuum.getDownQuarkEffectiveMass(), 
                                  vacuum.getStrangeQuarkEffectiveMass());

    cout << "inMediumSolution=" << inMedium.testSolution(1E-8) << "\n";
    cout << "Mu=" << inMedium.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << inMedium.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << inMedium.getStrangeQuarkEffectiveMass() << "GeV" << "\n";

    
    cout << "\n\n";

    double effMassU = inMedium.getUpQuarkEffectiveMass();
    double effMassD = inMedium.getDownQuarkEffectiveMass();
    double effMassS = inMedium.getStrangeQuarkEffectiveMass();

/*
    scatteringProcess process = UUUU;
    evaluateCrossSectionProcess12To34ToFile(parameters, T, 
                                            effChemPotU, effChemPotD, effChemPotS, 
                                            effMassU, effMassD, effMassS, 
                                            1E-8, process,  
                                            false, 1E-4,
                                            20);
*/
/*
    evaluateCrossSectionsKlevanskyPaper(parameters, T, 
                                        effChemPotU, effChemPotD, effChemPotS, 
                                        effMassU, effMassD, effMassS, 
                                        1E-8,
                                        false, 1E-4,
                                        200);
*/
/*
    double effMassU, effMassD, effMassS;

    double T = 0.250;
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSol = solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, T, 100, 1E-8, HYBRIDS);

    effMassU = finiteTempSol[int(finiteTempSol.size()-1)].getUpQuarkEffectiveMass();
    effMassD = finiteTempSol[int(finiteTempSol.size()-1)].getDownQuarkEffectiveMass();
    effMassS = finiteTempSol[int(finiteTempSol.size()-1)].getStrangeQuarkEffectiveMass();

    cout << "Mu=" << effMassU << "GeV" << "\t" 
         << "Md=" << effMassD << "GeV" << "\t" 
         << "Ms=" << effMassS << "GeV" << "\n";


    double chemPot = 0.100;
    vector<SU3NJL3DCutoffFixedChemPotTemp> inMediumSol = solveFromFiniteTemperatureToFiniteChemicalPotential(finiteTempSol[int(finiteTempSol.size()-1)], chemPot, 100, 1E-8, HYBRIDS);

    effMassU = inMediumSol[int(inMediumSol.size()-1)].getUpQuarkEffectiveMass();
    effMassD = inMediumSol[int(inMediumSol.size()-1)].getDownQuarkEffectiveMass();
    effMassS = inMediumSol[int(inMediumSol.size()-1)].getStrangeQuarkEffectiveMass();

    cout << "Mu=" << effMassU << "GeV" << "\t" 
         << "Md=" << effMassD << "GeV" << "\t" 
         << "Ms=" << effMassS << "GeV" << "\n";


    evaluateCrossSectionsPaperFiniteChemicalPotential(
        parameters, T, 
        0, 0, 0, 
        effMassU, effMassD, effMassS, 
        1E-8,
        false, 1E-4,
        200
    );

*/

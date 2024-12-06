#include <cmath>
#include <iostream>
#include <fstream>
#include "gsl_wrapper/root_solver_gsl.h"
#include "physics_utils/physical_constants.h"
#include "njl_model/line_integrals_3d_cutoff/OneFermionLineIntegral.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffBetaEqFixedTempRhoB.h"

using namespace std;



SU3NJL3DCutoffBetaEqFixedTempRhoB::SU3NJL3DCutoffBetaEqFixedTempRhoB(void* auxiliar)
{   
    parametersNJL = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->parametersNJL;
    electronMass = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->electronMass;
    temperature = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->temperature;
    baryonDensity = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->baryonDensity;

    upQuarkEffectiveMass = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->upQuarkEffectiveMass;
    downQuarkEffectiveMass = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->downQuarkEffectiveMass;
    strangeQuarkEffectiveMass = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->strangeQuarkEffectiveMass;
    upQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->upQuarkEffectiveChemicalPotential;
    downQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->downQuarkEffectiveChemicalPotential;
    strangeQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffBetaEqFixedTempRhoB *)(auxiliar))->strangeQuarkEffectiveChemicalPotential;
};



SU3NJL3DCutoffBetaEqFixedTempRhoB::SU3NJL3DCutoffBetaEqFixedTempRhoB(SU3NJL3DCutoffParameters parametersNJLAux, double electronMassAux, double temperatureAux, double baryonDensityAux)
{
    parametersNJL = parametersNJLAux;
    electronMass = electronMassAux;
    temperature = temperatureAux;
    baryonDensity = baryonDensityAux;
}



SU3NJL3DCutoffBetaEqFixedTempRhoB::SU3NJL3DCutoffBetaEqFixedTempRhoB(SU3NJL3DCutoffParameters parametersNJLAux, 
                                                                     double electronMassAux,
                                                                     double temperatureAux,
                                                                     double upQuarkEffectiveMassAux, 
                                                                     double downQuarkEffectiveMassAux, 
                                                                     double strangeQuarkEffectiveMassAux,
                                                                     double upQuarkEffectiveChemicalPotentialAux, 
                                                                     double downQuarkEffectiveChemicalPotentialAux, 
                                                                     double strangeQuarkEffectiveChemicalPotentialAux)
{
    parametersNJL = parametersNJLAux;
    electronMass = electronMassAux;
    temperature = temperatureAux;

    upQuarkEffectiveMass = upQuarkEffectiveMassAux;
    downQuarkEffectiveMass = downQuarkEffectiveMassAux;
    strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux;
    upQuarkEffectiveChemicalPotential = upQuarkEffectiveChemicalPotentialAux;
    downQuarkEffectiveChemicalPotential = downQuarkEffectiveChemicalPotentialAux;
    strangeQuarkEffectiveChemicalPotential = strangeQuarkEffectiveChemicalPotentialAux;
}



void SU3NJL3DCutoffBetaEqFixedTempRhoB::setSigmasAndDensitiesAndChemicalPotentials()
{
    //calculate the outputs of this solution
    NJLDimensionfulCouplings couplings = parametersNJL.getDimensionfulCouplings();

    NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();

    double Nc = parametersNJL.getNumberOfColours();
    
    double sigmaIntegralPrecision = parametersNJL.getSigmaIntegralPrecision();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

    double mE = getElectronMass();
    double T = getTemperature();
    
    
    double mU = getUpQuarkEffectiveMass();
    double mD = getDownQuarkEffectiveMass();
    double mS = getStrangeQuarkEffectiveMass();
    double effCPU = getUpQuarkEffectiveChemicalPotential();
    double effCPD = getDownQuarkEffectiveChemicalPotential();
    double effCPS = getStrangeQuarkEffectiveChemicalPotential();

    //Calculate the sigma field and density for each quark flavour
    double sigmaU = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCPU, mU, sigmaIntegralPrecision);
    double sigmaD = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCPD, mD, sigmaIntegralPrecision);
    double sigmaS = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCPS, mS, sigmaIntegralPrecision);

    double rhoU = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCPU, mU, thermoIntegralPrecision);
    double rhoD = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCPD, mD, thermoIntegralPrecision);
    double rhoS = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCPS, mS, thermoIntegralPrecision);

    //calculate quark chemical potential from the effective ones
    double cPU = SU3NJLQuarkChemicalPotential(couplings, effCPU, rhoU, rhoD, rhoS, sigmaU, sigmaD, sigmaS);
    double cPD = SU3NJLQuarkChemicalPotential(couplings, effCPD, rhoD, rhoS, rhoU, sigmaD, sigmaS, sigmaU);
    double cPS = SU3NJLQuarkChemicalPotential(couplings, effCPS, rhoS, rhoU, rhoD, sigmaS, sigmaU, sigmaD);

    //consider the electron in beta-equilibrium
    double cPE = cPD - cPU;
    double rhoE = fermionParticleDensity3DCutoff(reguScheme, cutoff, T, cPE, mE, thermoIntegralPrecision);

    setUpQuarkSigma(sigmaU);
    setDownQuarkSigma(sigmaD);
    setStrangeQuarkSigma(sigmaS);

    setUpQuarkDensity(rhoU);
    setDownQuarkDensity(rhoD);
    setStrangeQuarkDensity(rhoS);

    setUpQuarkChemicalPotential(cPU);
    setDownQuarkChemicalPotential(cPD);
    setStrangeQuarkChemicalPotential(cPS);
    setBaryonChemicalPotential();

    setElectronChemicalPotential(cPE);
    setElectronDensity(rhoE);
}



void SU3NJL3DCutoffBetaEqFixedTempRhoB::solve(double precision, 
                                              MultiRootFindingMethod method,
                                              double upQuarkEffectiveMassGuess, 
                                              double downQuarkEffectiveMassGuess, 
                                              double strangeQuarkEffectiveMassGuess, 
                                              double upQuarkEffectiveChemicalPotentialGuess, 
                                              double downQuarkEffectiveChemicalPotentialGuess, 
                                              double strangeQuarkEffectiveChemicalPotentialGuess)
{   
    //solve the model given the guesses
    double x[6];
    x[0] = upQuarkEffectiveMassGuess; 
    x[1] = downQuarkEffectiveMassGuess; 
    x[2] = strangeQuarkEffectiveMassGuess;
    x[3] = upQuarkEffectiveChemicalPotentialGuess;
    x[4] = downQuarkEffectiveChemicalPotentialGuess;
    x[5] = strangeQuarkEffectiveChemicalPotentialGuess;
    
    multiDimensionalRootFind(6, precision, &x[0], this, &SU3NJL3DCutoffGapEquationsBetaEquilibriumFixedTemperature, method);

    setUpQuarkEffectiveMass(x[0]);
    setDownQuarkEffectiveMass(x[1]);
    setStrangeQuarkEffectiveMass(x[2]);
    setUpQuarkEffectiveChemicalPotential(x[3]);
    setDownQuarkEffectiveChemicalPotential(x[4]);
    setStrangeQuarkEffectiveChemicalPotential(x[5]);

    //calculate outputs
    setSigmasAndDensitiesAndChemicalPotentials();
}



int SU3NJL3DCutoffGapEquationsBetaEquilibriumFixedTemperature(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{
	//define variables
    double mU = gsl_vector_get(x,0);
    double mD = gsl_vector_get(x,1);
    double mS = gsl_vector_get(x,2);
    double effCPU = gsl_vector_get(x,3);
    double effCPD = gsl_vector_get(x,4);
    double effCPS = gsl_vector_get(x,5);

    //read parameters from source object
    SU3NJL3DCutoffBetaEqFixedTempRhoB source(auxiliar);

    NJLDimensionfulCouplings couplings = source.getParametersNJL().getDimensionfulCouplings();

    double m0U = source.getParametersNJL().getUpQuarkCurrentMass();
    double m0D = source.getParametersNJL().getDownQuarkCurrentMass();
    double m0S = source.getParametersNJL().getStrangeQuarkCurrentMass();

    double mE = source.getElectronMass();
    double T = source.getTemperature();
    double rhoB = source.getBaryonDensity();


    //create an instance of the object with the guesses and calculate the outputs
    SU3NJL3DCutoffBetaEqFixedTempRhoB guess(source.getParametersNJL(), mE, T, mU, mD, mS, effCPU, effCPD, effCPS);
    guess.setSigmasAndDensitiesAndChemicalPotentials();

    double sigmaU = guess.getUpQuarkSigma();
    double sigmaD = guess.getDownQuarkSigma();
    double sigmaS = guess.getStrangeQuarkSigma();

    double rhoU = guess.getUpQuarkDensity();
    double rhoD = guess.getDownQuarkDensity();
    double rhoS = guess.getStrangeQuarkDensity();
    double rhoE = guess.getElectronDensity();

    double cPD = guess.getDownQuarkChemicalPotential();
    double cPS = guess.getStrangeQuarkChemicalPotential();
    

	//system of equations
    double f0 = SU3NJLNulledGapEquation(couplings, mU-m0U, sigmaU, sigmaD, sigmaS, rhoU, rhoD, rhoS);
    double f1 = SU3NJLNulledGapEquation(couplings, mD-m0D, sigmaD, sigmaS, sigmaU, rhoD, rhoS, rhoU);
    double f2 = SU3NJLNulledGapEquation(couplings, mS-m0S, sigmaS, sigmaU, sigmaD, rhoS, rhoU, rhoD);
    double f3 = rhoB - SU3BaryonDensity(rhoU, rhoD, rhoS);
	double f4 = cPD - cPS;
	double f5 = (2.0/3.0)*rhoU - (1./3)*rhoD - (1./3)*rhoS - rhoE;

	gsl_vector_set (f,0,f0);
	gsl_vector_set (f,1,f1);
	gsl_vector_set (f,2,f2);
	gsl_vector_set (f,3,f3);
	gsl_vector_set (f,4,f4);
	gsl_vector_set (f,5,f5);

	return GSL_SUCCESS;
}



bool SU3NJL3DCutoffBetaEqFixedTempRhoB::testSolution(double precision)
{   
    double x[6];
    x[0] = getUpQuarkEffectiveMass();
    x[1] = getDownQuarkEffectiveMass();
    x[2] = getStrangeQuarkEffectiveMass();
    x[3] = getUpQuarkEffectiveChemicalPotential();
    x[4] = getDownQuarkEffectiveChemicalPotential();
    x[5] = getStrangeQuarkEffectiveChemicalPotential();

    //the test below (gsl) resturns 0 if the sum_i abs(residual_i) < precision
    int gslTest = multiDimensionalRootFindTestResidual(6, precision, &x[0], this, &SU3NJL3DCutoffGapEquationsBetaEquilibriumFixedTemperature);

    if (gslTest==0){ return true; }
    else{ return false; }
}



double SU3NJL3DCutoffBetaEqFixedTempRhoB::calculatePressureQuarks(double vacuumPressureQuarks)
{   
    double T = getTemperature();
    double effMassU = getUpQuarkEffectiveMass();
    double effMassD = getDownQuarkEffectiveMass();
    double effMassS = getStrangeQuarkEffectiveMass();
    double effChemPotU = getUpQuarkEffectiveChemicalPotential();
    double effChemPotD = getDownQuarkEffectiveChemicalPotential();
    double effChemPotS = getStrangeQuarkEffectiveChemicalPotential();

    double pressureQuarks = SU3NJL3DCutoffPressure(parametersNJL, T, effMassU, effMassD, effMassS, effChemPotU, effChemPotD, effChemPotS);

    double pressure = pressureQuarks - vacuumPressureQuarks;

    return pressure;
}



double SU3NJL3DCutoffBetaEqFixedTempRhoB::calculatePressureElectrons(double vacuumPressureElectrons)
{   
    NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

    double massE = getElectronMass();
    double T = getTemperature();
    double chemPotE = getElectronChemicalPotential();

    double pressureElectrons = fermionPressure3DCutoff(reguScheme, cutoff, T, massE, chemPotE, thermoIntegralPrecision);

    double pressure = pressureElectrons - vacuumPressureElectrons;

    return pressure;
}



double SU3NJL3DCutoffBetaEqFixedTempRhoB::calculatePressure(double vacuumPressureQuarks, double vacuumPressureElectrons)
{
    double pressureQuarks = calculatePressureQuarks(vacuumPressureQuarks);
    double pressureElectrons = calculatePressureElectrons(vacuumPressureElectrons);

    double pressure = pressureQuarks + pressureElectrons;

    return pressure;
}



void SU3NJL3DCutoffBetaEqFixedTempRhoB::setBetaEqPressure(double vacuumPressureQuarks, double vacuumPressureElectrons)
{
    betaEqPressure = calculatePressure(vacuumPressureQuarks, vacuumPressureElectrons);
}



double SU3NJL3DCutoffBetaEqFixedTempRhoB::calculateEnergyDensityQuarks(double vacuumPressureQuarks)
{   
    double T = getTemperature();
    double effMassU = getUpQuarkEffectiveMass();
    double effMassD = getDownQuarkEffectiveMass();
    double effMassS = getStrangeQuarkEffectiveMass();
    double effChemPotU = getUpQuarkEffectiveChemicalPotential();
    double effChemPotD = getDownQuarkEffectiveChemicalPotential();
    double effChemPotS = getStrangeQuarkEffectiveChemicalPotential();

    double energyQuarks = SU3NJL3DCutoffEnergyDensity(parametersNJL, T, effMassU, effMassD, effMassS, effChemPotU, effChemPotD, effChemPotS);

    double energy = energyQuarks - (-vacuumPressureQuarks);

    return energy;
}



double SU3NJL3DCutoffBetaEqFixedTempRhoB::calculateEnergyDensityElectrons(double vacuumPressureElectrons)
{   
    NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

    double massE = getElectronMass();
    double T = getTemperature();
    double chemPotE = getElectronChemicalPotential();

    double energyElectrons = fermionEnergyDensity3DCutoff(reguScheme, cutoff, T, massE, chemPotE, chemPotE, thermoIntegralPrecision);

    double energy = energyElectrons - (-vacuumPressureElectrons);

    return energy;
}



double SU3NJL3DCutoffBetaEqFixedTempRhoB::calculateEnergyDensity(double vacuumPressureQuarks, double vacuumPressureElectrons)
{
    double energyQuarks = calculateEnergyDensityQuarks(vacuumPressureQuarks);

    double energyElectrons = calculateEnergyDensityElectrons(vacuumPressureElectrons);

    double energy = energyQuarks + energyElectrons;

    return energy;
}



void SU3NJL3DCutoffBetaEqFixedTempRhoB::setBetaEqEnergyDensity(double vacuumPressureQuarks, double vacuumPressureElectrons)
{
    betaEqEnergyDensity = calculateEnergyDensity(vacuumPressureQuarks, vacuumPressureElectrons);
}



double SU3NJL3DCutoffBetaEqFixedTempRhoB::calculateEntropyDensityQuarks()
{
    double T = getTemperature();
    double effMassU = getUpQuarkEffectiveMass();
    double effMassD = getDownQuarkEffectiveMass();
    double effMassS = getStrangeQuarkEffectiveMass();
    double effChemPotU = getUpQuarkEffectiveChemicalPotential();
    double effChemPotD = getDownQuarkEffectiveChemicalPotential();
    double effChemPotS = getStrangeQuarkEffectiveChemicalPotential();

    double entropyQuarks = SU3NJL3DCutoffEntropyDensity(parametersNJL, T, effMassU, effMassD, effMassS, effChemPotU, effChemPotD, effChemPotS);

    return entropyQuarks;
}



double SU3NJL3DCutoffBetaEqFixedTempRhoB::calculateEntropyDensityElectrons()
{
    NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

    double massE = getElectronMass();
    double T = getTemperature();
    double chemPotE = getElectronChemicalPotential();

    double entropyElectrons = fermionEntropyDensity3DCutoff(reguScheme, cutoff, T, massE, chemPotE, thermoIntegralPrecision);

    return entropyElectrons;
}



double SU3NJL3DCutoffBetaEqFixedTempRhoB::calculateEntropyDensity()
{
    double entropyQuarks = calculateEntropyDensityQuarks();

    double entropyElectrons = calculateEntropyDensityElectrons();

    double entropy = entropyQuarks + entropyElectrons;

    return entropy;
}



void SU3NJL3DCutoffBetaEqFixedTempRhoB::setBetaEqEntropyDensity()
{
    betaEqEntropyDensity = calculateEntropyDensity();
}



void SU3NJL3DCutoffBetaEqFixedTempRhoB::setBetaEqThermodynamics(double vacuumPressureQuarks, double vacuumPressureElectrons)
{
    setBetaEqPressure(vacuumPressureQuarks, vacuumPressureElectrons);
    setBetaEqEnergyDensity(vacuumPressureQuarks, vacuumPressureElectrons);
    setBetaEqEntropyDensity();
}


void writeSolutionsToFile(vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> solutions, string fileName, bool columnsDescription)
{
    int dataPrecision = 15;
    int colW = 25;

    std::ofstream fileSol;
    fileSol.open(fileName, std::fstream::in | std::ofstream::out | std::ios::trunc);
    fileSol.precision(dataPrecision);

    if ( columnsDescription )
    {
        fileSol.width(colW); fileSol << "rhoB[fm-3]";
        fileSol.width(colW); fileSol << "rhoU[fm-3]";
        fileSol.width(colW); fileSol << "rhoD[fm-3]";
        fileSol.width(colW); fileSol << "rhoS[fm-3]";
        fileSol.width(colW); fileSol << "rhoE[fm-3]";
        fileSol.width(colW); fileSol << "massU[MeV]";
        fileSol.width(colW); fileSol << "massD[MeV]";
        fileSol.width(colW); fileSol << "massS[MeV]";
        fileSol.width(colW); fileSol << "effChemPotU[MeV]";
        fileSol.width(colW); fileSol << "effChemPotD[MeV]";
        fileSol.width(colW); fileSol << "effChemPotS[MeV]";
        fileSol.width(colW); fileSol << "chemPotU[MeV]";
        fileSol.width(colW); fileSol << "chemPotD[MeV]";
        fileSol.width(colW); fileSol << "chemPotS[MeV]";
        fileSol.width(colW); fileSol << "baryonChemPot[MeV]";
        fileSol.width(colW); fileSol << "chemPotE[MeV]";
        fileSol.width(colW); fileSol << "pressure[MeVfm-3]";
        fileSol.width(colW); fileSol << "energyDensity[MeVfm-3]";
        fileSol.width(colW); fileSol << "entropyDensity[MeVfm-3]";
        fileSol << "\n";
    }

    for (int i = 0; i < int(solutions.size()); ++i)
    {
        fileSol.width(25);   fileSol << solutions[i].getBaryonDensity()/pow(hc_GeVfm,3);
        fileSol.width(25);   fileSol << solutions[i].getUpQuarkDensity()/pow(hc_GeVfm,3);
        fileSol.width(25);   fileSol << solutions[i].getDownQuarkDensity()/pow(hc_GeVfm,3);
        fileSol.width(25);   fileSol << solutions[i].getStrangeQuarkDensity()/pow(hc_GeVfm,3);
        fileSol.width(25);   fileSol << solutions[i].getElectronDensity()/pow(hc_GeVfm,3);

        fileSol.width(25);   fileSol << solutions[i].getUpQuarkEffectiveMass()*1000;
        fileSol.width(25);   fileSol << solutions[i].getDownQuarkEffectiveMass()*1000;
        fileSol.width(25);   fileSol << solutions[i].getStrangeQuarkEffectiveMass()*1000;

        fileSol.width(25);   fileSol << solutions[i].getUpQuarkEffectiveChemicalPotential()*1000;
        fileSol.width(25);   fileSol << solutions[i].getDownQuarkEffectiveChemicalPotential()*1000;
        fileSol.width(25);   fileSol << solutions[i].getStrangeQuarkEffectiveChemicalPotential()*1000;

        fileSol.width(25);   fileSol << solutions[i].getUpQuarkChemicalPotential()*1000;
        fileSol.width(25);   fileSol << solutions[i].getDownQuarkChemicalPotential()*1000;
        fileSol.width(25);   fileSol << solutions[i].getStrangeQuarkChemicalPotential()*1000;

        fileSol.width(25);   fileSol << solutions[i].getBaryonChemicalPotential()*1000;
        fileSol.width(25);   fileSol << solutions[i].getElectronChemicalPotential()*1000;

        fileSol.width(25);   fileSol << solutions[i].getBetaEqPressure()*(1000/pow(hc_GeVfm,3));
        fileSol.width(25);   fileSol << solutions[i].getBetaEqEnergyDensity()*(1000/pow(hc_GeVfm,3));
        fileSol.width(25);   fileSol << solutions[i].getBetaEqEntropyDensity()*(1000/pow(hc_GeVfm,3));

        fileSol << "\n";
    }

    fileSol.close();
}



void writeEOSToFile(vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> solutions, string fileName, bool columnsDescription)
{
    int dataPrecision = 15;
    int colW = 25;

    std::ofstream fileSol;
    fileSol.open(fileName, std::fstream::in | std::ofstream::out | std::ios::trunc);
    fileSol.precision(dataPrecision);

    if ( columnsDescription )
    {
        fileSol.width(colW); fileSol << "rhoB[fm-3]";
        fileSol.width(colW); fileSol << "energyDensity[MeVfm-3]";
        fileSol.width(colW); fileSol << "pressure[MeVfm-3]";
        fileSol.width(colW); fileSol << "baryonChemPot[MeV]";
        fileSol << "\n";
    }

    for (int i = 0; i < int(solutions.size()); ++i)
    {
        fileSol.width(25);   fileSol << solutions[i].getBaryonDensity()/pow(hc_GeVfm,3);
        fileSol.width(25);   fileSol << solutions[i].getBetaEqEnergyDensity()*(1000/pow(hc_GeVfm,3));
        fileSol.width(25);   fileSol << solutions[i].getBetaEqPressure()*(1000/pow(hc_GeVfm,3));
        fileSol.width(25);   fileSol << solutions[i].getBaryonChemicalPotential()*1000;
        fileSol << "\n";
    }

    fileSol.close();
}


void writeEOSToFile(vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> solutions, string fileName, bool columnsDescription, double minimumBaryonDensity)
{
    int dataPrecision = 15;
    int colW = 25;

    std::ofstream fileSol;
    fileSol.open(fileName, std::fstream::in | std::ofstream::out | std::ios::trunc);
    fileSol.precision(dataPrecision);

    if ( columnsDescription )
    {
        fileSol.width(colW); fileSol << "rhoB[fm-3]";
        fileSol.width(colW); fileSol << "energyDensity[MeVfm-3]";
        fileSol.width(colW); fileSol << "pressure[MeVfm-3]";
        fileSol.width(colW); fileSol << "baryonChemPot[MeV]";
        fileSol << "\n";
    }

    for (int i = 0; i < int(solutions.size()); ++i)
    {   
        if ( solutions[i].getBaryonDensity()>minimumBaryonDensity )
        {
            fileSol.width(25);   fileSol << solutions[i].getBaryonDensity()/pow(hc_GeVfm,3);
            fileSol.width(25);   fileSol << solutions[i].getBetaEqEnergyDensity()*(1000/pow(hc_GeVfm,3));
            fileSol.width(25);   fileSol << solutions[i].getBetaEqPressure()*(1000/pow(hc_GeVfm,3));
            fileSol.width(25);   fileSol << solutions[i].getBaryonChemicalPotential()*1000;
            fileSol << "\n";
        }
    }

    fileSol.close();
}


int SU3NJL3DCutoffChiralTransitionPointBetaEquilibriumFixedTemperature(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{
    //define variables
    double mU_broken = gsl_vector_get(x,0);
    double mD_broken = gsl_vector_get(x,1);
    double mS_broken = gsl_vector_get(x,2);
    double effCPU_broken = gsl_vector_get(x,3);
    double effCPD_broken = gsl_vector_get(x,4);
    double effCPS_broken = gsl_vector_get(x,5);

    double mU_restored = gsl_vector_get(x,6);
    double mD_restored = gsl_vector_get(x,7);
    double mS_restored = gsl_vector_get(x,8);
    double effCPU_restored = gsl_vector_get(x,9);
    double effCPD_restored = gsl_vector_get(x,10);
    double effCPS_restored = gsl_vector_get(x,11);


    //read parameters from source object
    SU3NJL3DCutoffBetaEqFixedTempRhoB source(auxiliar);

    NJLDimensionfulCouplings couplings = source.getParametersNJL().getDimensionfulCouplings();
    double m0U = source.getParametersNJL().getUpQuarkCurrentMass();
    double m0D = source.getParametersNJL().getDownQuarkCurrentMass();
    double m0S = source.getParametersNJL().getStrangeQuarkCurrentMass();
    double mE = source.getElectronMass();
    double T = source.getTemperature();


    //Broken phase
    SU3NJL3DCutoffBetaEqFixedTempRhoB brokenPhase(source.getParametersNJL(), 
                                                  mE, 
                                                  T,
                                                  mU_broken, 
                                                  mD_broken, 
                                                  mS_broken, 
                                                  effCPU_broken, 
                                                  effCPD_broken, 
                                                  effCPS_broken);
    brokenPhase.setSigmasAndDensitiesAndChemicalPotentials();

    double sigmaU_broken = brokenPhase.getUpQuarkSigma();
    double sigmaD_broken = brokenPhase.getDownQuarkSigma();
    double sigmaS_broken = brokenPhase.getStrangeQuarkSigma();

    //get all the necessary quantities
    double rhoU_broken = brokenPhase.getUpQuarkDensity();
    double rhoD_broken = brokenPhase.getDownQuarkDensity();
    double rhoS_broken = brokenPhase.getStrangeQuarkDensity();
    double rhoE_broken = brokenPhase.getElectronDensity();

    double cPD_broken = brokenPhase.getDownQuarkChemicalPotential();
    double cPS_broken = brokenPhase.getStrangeQuarkChemicalPotential();
    double baryonCP_broken = brokenPhase.getBaryonChemicalPotential();

    double pressure_broken = brokenPhase.calculatePressure(0.0, 0.0);//the vacuum values can be any value


    //Restored phase
    SU3NJL3DCutoffBetaEqFixedTempRhoB restoredPhase(source.getParametersNJL(), 
                                                    mE, 
                                                    T,
                                                    mU_restored, 
                                                    mD_restored, 
                                                    mS_restored, 
                                                    effCPU_restored, 
                                                    effCPD_restored, 
                                                    effCPS_restored);
    restoredPhase.setSigmasAndDensitiesAndChemicalPotentials();

    double sigmaU_restored = restoredPhase.getUpQuarkSigma();
    double sigmaD_restored = restoredPhase.getDownQuarkSigma();
    double sigmaS_restored = restoredPhase.getStrangeQuarkSigma();

    //get all the necessary quantities
    double rhoU_restored = restoredPhase.getUpQuarkDensity();
    double rhoD_restored = restoredPhase.getDownQuarkDensity();
    double rhoS_restored = restoredPhase.getStrangeQuarkDensity();
    double rhoE_restored = restoredPhase.getElectronDensity();

    double cPD_restored = restoredPhase.getDownQuarkChemicalPotential();
    double cPS_restored = restoredPhase.getStrangeQuarkChemicalPotential();
    double baryonCP_restored = restoredPhase.getBaryonChemicalPotential();

    double pressure_restored = restoredPhase.calculatePressure(0.0, 0.0);//the vacuum values can be any value


    //equations in the broken phase
    double f0 = SU3NJLNulledGapEquation(couplings, mU_broken-m0U, sigmaU_broken, sigmaD_broken, sigmaS_broken, rhoU_broken, rhoD_broken, rhoS_broken);
    double f1 = SU3NJLNulledGapEquation(couplings, mD_broken-m0D, sigmaD_broken, sigmaS_broken, sigmaU_broken, rhoD_broken, rhoS_broken, rhoU_broken);
    double f2 = SU3NJLNulledGapEquation(couplings, mS_broken-m0S, sigmaS_broken, sigmaU_broken, sigmaD_broken, rhoS_broken, rhoU_broken, rhoD_broken);
    double f3 = cPD_broken - cPS_broken;
    double f4 = (2.0/3.0)*rhoU_broken - (1.0/3.0)*rhoD_broken - (1.0/3.0)*rhoS_broken - rhoE_broken;

    //equations in the restored phase
    double f5 = SU3NJLNulledGapEquation(couplings, mU_restored-m0U, sigmaU_restored, sigmaD_restored, sigmaS_restored, rhoU_restored, rhoD_restored, rhoS_restored);
    double f6 = SU3NJLNulledGapEquation(couplings, mD_restored-m0D, sigmaD_restored, sigmaS_restored, sigmaU_restored, rhoD_restored, rhoS_restored, rhoU_restored);
    double f7 = SU3NJLNulledGapEquation(couplings, mS_restored-m0S, sigmaS_restored, sigmaU_restored, sigmaD_restored, rhoS_restored, rhoU_restored, rhoD_restored);
    double f8 = cPD_restored - cPS_restored;
    double f9 = (2.0/3.0)*rhoU_restored - (1.0/3.0)*rhoD_restored - (1.0/3.0)*rhoS_restored - rhoE_restored;
    
    //Maxwell construction
    double f10 = baryonCP_broken - baryonCP_restored;
    double f11 = pressure_broken - pressure_restored;


    gsl_vector_set (f,0,f0);
    gsl_vector_set (f,1,f1);
    gsl_vector_set (f,2,f2);
    gsl_vector_set (f,3,f3);
    gsl_vector_set (f,4,f4);
    gsl_vector_set (f,5,f5);
    gsl_vector_set (f,6,f6);
    gsl_vector_set (f,7,f7);
    gsl_vector_set (f,8,f8);
    gsl_vector_set (f,9,f9);
    gsl_vector_set (f,10,f10);
    gsl_vector_set (f,11,f11);

    return GSL_SUCCESS;
}



vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> findChiralTransitionPointsFixedTemperature(vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> solutions, double precision, MultiRootFindingMethod method)
{
    //run check that solutions have the same temperature
    for (int i = 0; i < int(solutions.size()); ++i)
    {
        if ( fabs(solutions[0].getTemperature()-solutions[i].getTemperature())>0 ){ cout << "Temperature is not fixed in this vector of solutions!\n"; abort(); }
    }


    //finding guesses for transition point
    int first = 0;
    for (int i = 1; i < int(solutions.size()); ++i)
    {
        double delta = solutions[i-1].getBaryonChemicalPotential() - solutions[i].getBaryonChemicalPotential();
        if ( delta>0 )
        {
            first = i;
            //cout << first << "\n";
            break;
        }
    }
    int second = 0;
    if ( first>0 )
    {
        for (int i = first; i < int(solutions.size()); ++i)
        {
            double delta = solutions[i-1].getBaryonChemicalPotential() - solutions[i].getBaryonChemicalPotential();
            if ( delta<0 )
            {
                second = i;
                //cout << second << "\n";
                break;
            }
        }
    }


/*
    //other guessing strategy
    int guess_broken = 0;
    bool run_broken = true;
    for (int i = 0; i < first; ++i)
    {   
        for (int j = second; j < int(solutions.size()); ++j)
        {
            if ( solutions[i].getBaryonChemicalPotential()>solutions[j].getBaryonChemicalPotential() )
            {   
                guess_broken = i;
                run_broken = false;
                break;
            }
        }
        if (run_broken==false){ break; }
    }
    int guess_restored = 0;
    bool run_restored = true;
    for (int i = second; i < int(solutions.size()); ++i)
    {   
        for (int j = 0; j < first; ++j)
        {
            if ( solutions[i].getBetaEqPressure()>solutions[j].getBetaEqPressure() )
            {   
                guess_restored = i;
                run_restored = false;
                break;
            }
        }
        if (run_restored==false){ break; }
    }
*/


    vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> transitionPoints;
    if ( first>0 && second>0 )
    {

    //estimate guesses
    double mU_broken = solutions[first].getUpQuarkEffectiveMass();
    double mD_broken = solutions[first].getDownQuarkEffectiveMass();
    double mS_broken = solutions[first].getStrangeQuarkEffectiveMass();
    double effCPU_broken = solutions[first].getUpQuarkEffectiveChemicalPotential();
    double effCPD_broken = solutions[first].getDownQuarkEffectiveChemicalPotential();
    double effCPS_broken = solutions[first].getStrangeQuarkEffectiveChemicalPotential();
    double mU_restored = solutions[second].getUpQuarkEffectiveMass();
    double mD_restored = solutions[second].getDownQuarkEffectiveMass();
    double mS_restored = solutions[second].getStrangeQuarkEffectiveMass();
    double effCPU_restored = solutions[second].getUpQuarkEffectiveChemicalPotential();
    double effCPD_restored = solutions[second].getDownQuarkEffectiveChemicalPotential();
    double effCPS_restored = solutions[second].getStrangeQuarkEffectiveChemicalPotential();


    //solve system
    double x[12];
    x[0] = mU_broken; 
    x[1] = mD_broken;  
    x[2] = mS_broken;  
    x[3] = 0.5*effCPU_broken + 0.5*effCPU_restored;
    x[4] = 0.5*effCPD_broken + 0.5*effCPD_restored;
    x[5] = 0.5*effCPS_broken + 0.5*effCPS_restored;
    x[6] = mU_restored; 
    x[7] = mD_restored;  
    x[8] = mS_restored;  
    x[9] = 0.5*effCPU_broken + 0.5*effCPU_restored;
    x[10] = 0.5*effCPD_broken + 0.5*effCPD_restored;
    x[11] = 0.5*effCPS_broken + 0.5*effCPS_restored;
    
    SU3NJL3DCutoffBetaEqFixedTempRhoB aux(solutions[0].getParametersNJL(), solutions[0].getElectronMass(), solutions[0].getTemperature(), 0.0);
    multiDimensionalRootFind(12, precision, &x[0], &aux, &SU3NJL3DCutoffChiralTransitionPointBetaEquilibriumFixedTemperature, method);

    mU_broken = x[0]; 
    mD_broken = x[1];
    mS_broken = x[2]; 
    effCPU_broken = x[3];
    effCPD_broken = x[4];
    effCPS_broken = x[5];
    mU_restored = x[6]; 
    mD_restored = x[7];  
    mS_restored = x[8];  
    effCPU_restored = x[9]; 
    effCPD_restored = x[10]; 
    effCPS_restored = x[11];

/*
    cout << "Broken phase:\n";
    cout << mU_broken << "\t" << mD_broken << "\t" << mS_broken << "\n";
    cout << effCPU_broken << "\t" << effCPD_broken << "\t" << effCPS_broken << "\n";
    cout << "Restored phase:\n";
    cout << mU_restored << "\t" << mD_restored << "\t" << mS_restored << "\n";
    cout << effCPU_restored << "\t" << effCPD_restored << "\t" << effCPS_restored << "\n";
*/

    SU3NJL3DCutoffBetaEqFixedTempRhoB broken(solutions[0].getParametersNJL(), 
                                             solutions[0].getElectronMass(), 
                                             solutions[0].getTemperature(), 
                                             mU_broken, 
                                             mD_broken, 
                                             mS_broken, 
                                             effCPU_broken, 
                                             effCPD_broken, 
                                             effCPS_broken);   
    broken.setSigmasAndDensitiesAndChemicalPotentials();
    double baryonDensity_broken = SU3BaryonDensity(broken.getUpQuarkDensity(), broken.getDownQuarkDensity(), broken.getStrangeQuarkDensity());
    broken.setBaryonDensity(baryonDensity_broken);


    SU3NJL3DCutoffBetaEqFixedTempRhoB restored(solutions[0].getParametersNJL(), 
                                               solutions[0].getElectronMass(), 
                                               solutions[0].getTemperature(), 
                                               mU_restored, 
                                               mD_restored, 
                                               mS_restored, 
                                               effCPU_restored, 
                                               effCPD_restored, 
                                               effCPS_restored);
    restored.setSigmasAndDensitiesAndChemicalPotentials();
    double baryonDensity_restored = SU3BaryonDensity(restored.getUpQuarkDensity(), restored.getDownQuarkDensity(), restored.getStrangeQuarkDensity());
    restored.setBaryonDensity(baryonDensity_restored);


    transitionPoints.push_back(broken);
    transitionPoints.push_back(restored);

    }


    return transitionPoints;
}



void addVacuumSolution(SU3NJL3DCutoffVacuum vacuum, double electronMass, double pressureVac, double pressureVacElec, vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> &solutions)
{
    SU3NJL3DCutoffBetaEqFixedTempRhoB aux(vacuum.getParametersNJL(), electronMass, 0.0,
                                          vacuum.getUpQuarkEffectiveMass(),
                                          vacuum.getDownQuarkEffectiveMass(),
                                          vacuum.getStrangeQuarkEffectiveMass(), 0.0, 0.0, 0.0);
    aux.setBaryonDensity(0.0);
    aux.setSigmasAndDensitiesAndChemicalPotentials();
    aux.setBetaEqThermodynamics(pressureVac, pressureVacElec);

    solutions.push_back(aux);
}


vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> 
calculateZeroTemperatureSolutions(SU3NJL3DCutoffVacuum vacuum, 
                                  double minimumBaryonDensity, double maximumBaryonDensity, int numberOfPoints, 
                                  double gapPrecision, MultiRootFindingMethod method)
{
    double pressureVacuum = vacuum.calculatePressure();
    double energyVacuum = vacuum.calculateEnergyDensity();
    cout << "pressureVacuum=" << pressureVacuum << "\n";
    cout << "energyDensityVacuum=" << energyVacuum << "\n";

    double pressureVacuumElectron = vacuum.calculateVacuumPressureElectrons(electronMass_GeV);
    cout << "pressureVacuumElectron=" << pressureVacuumElectron << "\n";


    vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> betaEqSolutions;
    addVacuumSolution(vacuum, electronMass_GeV, pressureVacuum, pressureVacuumElectron, betaEqSolutions);


    //guesses for i=0
    double mUGuess = vacuum.getUpQuarkEffectiveMass();
    double mDGuess = vacuum.getDownQuarkEffectiveMass();
    double mSGuess = vacuum.getStrangeQuarkEffectiveMass();
    double effCPUGuess = mUGuess + 1E-3;
    double effCPDGuess = mUGuess + 1E-3;
    double effCPSGuess = mUGuess + 1E-3;


    double temperature = 0.0;
    double delta = (maximumBaryonDensity-minimumBaryonDensity)/(numberOfPoints-1);
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        //step in density
        double rhoB = minimumBaryonDensity + i*delta;

        //create object
        SU3NJL3DCutoffBetaEqFixedTempRhoB betaEq(vacuum.getParametersNJL(), electronMass_GeV, temperature, rhoB);

        //find quark masses and effective chemical potential
        betaEq.solve(gapPrecision, method, mUGuess, mDGuess, mSGuess, effCPUGuess, effCPDGuess, effCPSGuess);

        //guesses for next step
        mUGuess = betaEq.getUpQuarkEffectiveMass();
        mDGuess = betaEq.getDownQuarkEffectiveMass();
        mSGuess = betaEq.getStrangeQuarkEffectiveMass();
        effCPUGuess = betaEq.getUpQuarkEffectiveChemicalPotential();
        effCPDGuess = betaEq.getDownQuarkEffectiveChemicalPotential();
        effCPSGuess = betaEq.getStrangeQuarkEffectiveChemicalPotential();

        //calculate thermodynamics
        betaEq.setBetaEqThermodynamics(pressureVacuum, pressureVacuumElectron);


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


    return betaEqSolutions;
}


void writeBetaEquilibriumEOSAtZeroTemperatureToFile(SU3NJL3DCutoffVacuum vacuum, 
                                                    double minimumBaryonDensity, double maximumBaryonDensity, int numberOfPoints, 
                                                    double gapPrecision, MultiRootFindingMethod method,
                                                    string fileName)
{
    vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> betaEqSolutions = 
    calculateZeroTemperatureSolutions(vacuum, minimumBaryonDensity, maximumBaryonDensity, numberOfPoints, gapPrecision, method);


    //save all information in file
    writeSolutionsToFile(betaEqSolutions, "solutions.dat", true);


    //find chiral transition if it exists
    vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> transitionPoints = findChiralTransitionPointsFixedTemperature(betaEqSolutions, gapPrecision, method);


    //save EOS to file: if it has first order phase transition, save only after restoration 
    if ( int(transitionPoints.size())>0 )
    {   
        double minRhoB = transitionPoints[1].getBaryonDensity();
        writeEOSToFile(betaEqSolutions, fileName, true, minRhoB);
    }
    else
    {
        writeEOSToFile(betaEqSolutions, fileName, true);
    }
}

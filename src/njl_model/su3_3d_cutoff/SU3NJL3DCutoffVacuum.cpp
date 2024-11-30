#include <cmath>
#include <iostream>
#include "gsl_wrapper/root_solver_gsl.h"
#include "OneFermionLineIntegral.h"
#include "SU3NJL3DCutoff.h"
#include "SU3NJL3DCutoffVacuum.h"

using namespace std;


SU3NJL3DCutoffVacuum::SU3NJL3DCutoffVacuum(SU3NJL3DCutoffParameters parametersNJLAux)
{
	parametersNJL = parametersNJLAux;
}


SU3NJL3DCutoffVacuum::SU3NJL3DCutoffVacuum(void* auxiliar)
{	
	parametersNJL = ((class SU3NJL3DCutoffVacuum *)(auxiliar))->parametersNJL;

	upQuarkEffectiveMass = ((class SU3NJL3DCutoffVacuum *)(auxiliar))->upQuarkEffectiveMass;
	downQuarkEffectiveMass = ((class SU3NJL3DCutoffVacuum *)(auxiliar))->downQuarkEffectiveMass;
	strangeQuarkEffectiveMass = ((class SU3NJL3DCutoffVacuum *)(auxiliar))->strangeQuarkEffectiveMass;
};


void SU3NJL3DCutoffVacuum::solve(double precision, MultiRootFindingMethod method, double upQuarkEffectiveMassGuess, double downQuarkEffectiveMassGuess, double strangeQuarkEffectiveMassGuess)
{	
	double x[3];
    x[0] = upQuarkEffectiveMassGuess; 
    x[1] = downQuarkEffectiveMassGuess; 
    x[2] = strangeQuarkEffectiveMassGuess; 
    
    multiDimensionalRootFind(3, precision, &x[0], this, &SU3NJLGapEquationsVacuum, method);

	setUpQuarkEffectiveMass(x[0]);
	setDownQuarkEffectiveMass(x[1]);
	setStrangeQuarkEffectiveMass(x[2]);
}


int SU3NJLGapEquationsVacuum(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{   
    //define variables
    double mU = gsl_vector_get(x,0);
    double mD = gsl_vector_get(x,1);
    double mS = gsl_vector_get(x,2);

    //define parameters
    SU3NJL3DCutoffVacuum solution(auxiliar);

    NJLDimensionfulCouplings couplings = solution.getParametersNJL().getDimensionfulCouplings();
    
    NJL3DCutoffRegularizationScheme reguScheme = solution.getParametersNJL().getNJL3DCutoffRegularizationScheme();
    double cutoff = solution.getParametersNJL().getThreeMomentumCutoff();

    double Nc = solution.getParametersNJL().getNumberOfColours();

    double sigmaIntegralPrecision = solution.getParametersNJL().getSigmaIntegralPrecision();

    double m0U = solution.getParametersNJL().getUpQuarkCurrentMass();
    double m0D = solution.getParametersNJL().getDownQuarkCurrentMass();
    double m0S = solution.getParametersNJL().getStrangeQuarkCurrentMass();


    //calculate quark condensates
    double T = 0.0;
    double effCP = 0.0;
    double sigmaU = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCP, mU, sigmaIntegralPrecision);
   	double sigmaD = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCP, mD, sigmaIntegralPrecision);
    double sigmaS = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCP, mS, sigmaIntegralPrecision);

    //gap equations
    double f0 = SU3NJLNulledGapEquation(couplings, mU-m0U, sigmaU, sigmaD, sigmaS, 0.0, 0.0, 0.0);
    double f1 = SU3NJLNulledGapEquation(couplings, mD-m0D, sigmaD, sigmaS, sigmaU, 0.0, 0.0, 0.0);
    double f2 = SU3NJLNulledGapEquation(couplings, mS-m0S, sigmaS, sigmaU, sigmaD, 0.0, 0.0, 0.0);

    gsl_vector_set (f, 0, f0);
    gsl_vector_set (f, 1, f1);
    gsl_vector_set (f, 2, f2);

    return GSL_SUCCESS;
}


bool SU3NJL3DCutoffVacuum::testSolution(double precision)
{   
    double x[3];
    x[0] = getUpQuarkEffectiveMass();
    x[1] = getDownQuarkEffectiveMass();
    x[2] = getStrangeQuarkEffectiveMass();
    
    //the test below (gsl) resturns 0 if the sum_i abs(residual_i) < precision
    int gslTest = multiDimensionalRootFindTestResidual(3, precision, &x[0], this, &SU3NJLGapEquationsVacuum);

    if (gslTest==0){ return true; }
    else{ return false; }
}


double SU3NJL3DCutoffVacuum::calculatePressure()
{
    double pressure = SU3NJL3DCutoffPressure(parametersNJL, 0.0, upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 0.0, 0.0, 0.0);
    return pressure;
}


double SU3NJL3DCutoffVacuum::calculateEnergyDensity()
{
    double energy = SU3NJL3DCutoffEnergyDensity(parametersNJL, 0.0, upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 0.0, 0.0, 0.0);
    return energy;
}


double SU3NJL3DCutoffVacuum::calculateVacuumPressureElectrons(double electronMass)
{   
    NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

    double vacuumPressureElectrons = fermionPressure3DCutoff(reguScheme, cutoff, 0.0, electronMass, 0.0, thermoIntegralPrecision);

    return vacuumPressureElectrons;
}


SU3NJL3DCutoffMeson SU3NJL3DCutoffVacuum::calculateMesonMassAndWidth(mesonState meson, double precision, MultiRootFindingMethod method, double mesonMassGuess, double mesonWidthGuess)
{   
    double temperature = 0.0;
    double effChemPot = 0.0;
    double mesonPropagatorPrecision = parametersNJL.getSigmaIntegralPrecision();

    SU3NJL3DCutoffMeson mesonAux(parametersNJL, temperature, 
                                 effChemPot, effChemPot, effChemPot, 
                                 upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                                 mesonPropagatorPrecision, meson);

    mesonAux.calculateMesonMassAndWidth(precision, method, mesonMassGuess, mesonWidthGuess);

    return mesonAux;
}




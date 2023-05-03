#include <cmath>
#include <iostream>
#include "rootSolverGSL.h"
#include "OneFermionLineIntegral.h"
#include "SU3NJL3DCutoff.h"
#include "SU3NJL3DCutoffFixedChemPotTemp.h"

using namespace std;


SU3NJL3DCutoffFixedChemPotTemp::SU3NJL3DCutoffFixedChemPotTemp(void* auxiliar)
{	
	parametersNJL = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->parametersNJL;

	temperature = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->temperature;

	upQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->upQuarkEffectiveChemicalPotential;
	downQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->downQuarkEffectiveChemicalPotential;
	strangeQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->strangeQuarkEffectiveChemicalPotential;

	upQuarkEffectiveMass = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->upQuarkEffectiveMass;
	downQuarkEffectiveMass = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->downQuarkEffectiveMass;
	strangeQuarkEffectiveMass = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->strangeQuarkEffectiveMass;
};


SU3NJL3DCutoffFixedChemPotTemp::SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffParameters parametersNJLAux, 
															   double temperatureAux, 
															   double upQuarkEffectiveChemicalPotentialAux,
															   double downQuarkEffectiveChemicalPotentialAux,
															   double strangeQuarkEffectiveChemicalPotentialAux)
{
	parametersNJL = parametersNJLAux;

	temperature = temperatureAux;

	upQuarkEffectiveChemicalPotential = upQuarkEffectiveChemicalPotentialAux;
	downQuarkEffectiveChemicalPotential = downQuarkEffectiveChemicalPotentialAux;
	strangeQuarkEffectiveChemicalPotential = strangeQuarkEffectiveChemicalPotentialAux;
}


void SU3NJL3DCutoffFixedChemPotTemp::solve(double precision, 
										   MultiRootFindingMethod method, 
                                           double upQuarkEffectiveMassGuess, 
                                           double downQuarkEffectiveMassGuess, 
                                           double strangeQuarkEffectiveMassGuess)
{	
	double x[3];
    x[0] = upQuarkEffectiveMassGuess; 
    x[1] = downQuarkEffectiveMassGuess; 
    x[2] = strangeQuarkEffectiveMassGuess;
    
    multiDimensionalRootFind(3, precision, &x[0], this, &SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature, method);

	setUpQuarkEffectiveMass(x[0]);
	setDownQuarkEffectiveMass(x[1]);
	setStrangeQuarkEffectiveMass(x[2]);
}


int SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{
	//define variables
    double mU = gsl_vector_get(x,0);
    double mD = gsl_vector_get(x,1);
    double mS = gsl_vector_get(x,2);


    //define parameters
    SU3NJL3DCutoffFixedChemPotTemp solution(auxiliar);

    NJLDimensionfulCouplings couplings = solution.getParametersNJL().getDimensionfulCouplings();

    NJL3DCutoffRegularizationScheme reguScheme = solution.getParametersNJL().getNJL3DCutoffRegularizationScheme();
    double cutoff = solution.getParametersNJL().getThreeMomentumCutoff();

    double Nc = solution.getParametersNJL().getNumberOfColours();

    double sigmaIntegralPrecision = solution.getParametersNJL().getSigmaIntegralPrecision();

    double m0U = solution.getParametersNJL().getUpQuarkCurrentMass();
    double m0D = solution.getParametersNJL().getDownQuarkCurrentMass();
    double m0S = solution.getParametersNJL().getStrangeQuarkCurrentMass();

    double T = solution.getTemperature();
    double effCPU = solution.getUpQuarkEffectiveChemicalPotential();
    double effCPD = solution.getDownQuarkEffectiveChemicalPotential();
    double effCPS = solution.getStrangeQuarkEffectiveChemicalPotential();


    //Calculate the sigma field and density for each quark flavour
    double sigmaU = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCPU, mU, sigmaIntegralPrecision);
    double sigmaD = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCPD, mD, sigmaIntegralPrecision);
    double sigmaS = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCPS, mS, sigmaIntegralPrecision);

	//system of equations
    double f0 = SU3NJLNulledGapEquation(couplings, mU-m0U, sigmaU, sigmaD, sigmaS, 0.0, 0.0, 0.0);
    double f1 = SU3NJLNulledGapEquation(couplings, mD-m0D, sigmaD, sigmaS, sigmaU, 0.0, 0.0, 0.0);
    double f2 = SU3NJLNulledGapEquation(couplings, mS-m0S, sigmaS, sigmaU, sigmaD, 0.0, 0.0, 0.0);
    
    /*
    double thermoIntegralPrecision = solution.getParametersNJL().getThermoIntegralPrecision();

    double rhoU = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCPU, mU, thermoIntegralPrecision);
    double rhoD = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCPD, mD, thermoIntegralPrecision);
    double rhoS = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCPS, mS, thermoIntegralPrecision);

   	//system of equations
    double f0 = SU3NJLNulledGapEquation(couplings, mU-m0U, sigmaU, sigmaD, sigmaS, rhoU, rhoD, rhoS);
    double f1 = SU3NJLNulledGapEquation(couplings, mD-m0D, sigmaD, sigmaS, sigmaU, rhoD, rhoS, rhoU);
    double f2 = SU3NJLNulledGapEquation(couplings, mS-m0S, sigmaS, sigmaU, sigmaD, rhoS, rhoU, rhoD);
	*/

	gsl_vector_set (f, 0, f0);
	gsl_vector_set (f, 1, f1);
	gsl_vector_set (f, 2, f2);

	return GSL_SUCCESS;
}


bool SU3NJL3DCutoffFixedChemPotTemp::testSolution(double precision)
{   
    double x[3];
    x[0] = getUpQuarkEffectiveMass();
    x[1] = getDownQuarkEffectiveMass();
    x[2] = getStrangeQuarkEffectiveMass();

    //the test below (gsl) resturns 0 if the sum_i abs(residual_i) < precision
    int gslTest = multiDimensionalRootFindTestResidual(3, precision, &x[0], this, &SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature);

    if (gslTest==0){ return true; }
    else{ return false; }
}



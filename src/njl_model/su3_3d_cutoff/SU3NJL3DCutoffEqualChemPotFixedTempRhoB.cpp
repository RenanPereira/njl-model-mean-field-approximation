#include <cmath>
#include <iostream>
#include "gsl_wrapper/root_solver_gsl.h"
#include "njl_model/line_integrals_3d_cutoff/OneFermionLineIntegral.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffEqualChemPotFixedTempRhoB.h"

using namespace std;


SU3NJL3DCutoffEqualChemPotFixedTempRhoB::SU3NJL3DCutoffEqualChemPotFixedTempRhoB(void* auxiliar)
{	
	parametersNJL = ((class SU3NJL3DCutoffEqualChemPotFixedTempRhoB *)(auxiliar))->parametersNJL;

	upQuarkEffectiveMass = ((class SU3NJL3DCutoffEqualChemPotFixedTempRhoB *)(auxiliar))->upQuarkEffectiveMass;
	downQuarkEffectiveMass = ((class SU3NJL3DCutoffEqualChemPotFixedTempRhoB *)(auxiliar))->downQuarkEffectiveMass;
	strangeQuarkEffectiveMass = ((class SU3NJL3DCutoffEqualChemPotFixedTempRhoB *)(auxiliar))->strangeQuarkEffectiveMass;
	quarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffEqualChemPotFixedTempRhoB *)(auxiliar))->quarkEffectiveChemicalPotential;

	temperature = ((class SU3NJL3DCutoffEqualChemPotFixedTempRhoB *)(auxiliar))->temperature;
	baryonDensity = ((class SU3NJL3DCutoffEqualChemPotFixedTempRhoB *)(auxiliar))->baryonDensity;
};


SU3NJL3DCutoffEqualChemPotFixedTempRhoB::SU3NJL3DCutoffEqualChemPotFixedTempRhoB(SU3NJL3DCutoffParameters parametersNJLAux)
{
	parametersNJL = parametersNJLAux;
}


void SU3NJL3DCutoffEqualChemPotFixedTempRhoB::solve(double precision, MultiRootFindingMethod method, 
                                                    double upQuarkEffectiveMassGuess, 
                                                    double downQuarkEffectiveMassGuess, 
                                                    double strangeQuarkEffectiveMassGuess, 
                                                    double effectiveChemicalPotentialGuess)
{	
	double x[4];
    x[0] = upQuarkEffectiveMassGuess; 
    x[1] = downQuarkEffectiveMassGuess; 
    x[2] = strangeQuarkEffectiveMassGuess;
    x[3] = effectiveChemicalPotentialGuess;
    
    multiDimensionalRootFind(4, precision, &x[0], this, &SU3NJL3DCutoffGapEquationsEqualChemicalPotentialFixedTemperatureBaryonDensity, method);

	setUpQuarkEffectiveMass(x[0]);
	setDownQuarkEffectiveMass(x[1]);
	setStrangeQuarkEffectiveMass(x[2]);
	setQuarkEffectiveChemicalPotential(x[3]);
}


int SU3NJL3DCutoffGapEquationsEqualChemicalPotentialFixedTemperatureBaryonDensity(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{
	//define variables
    double mU = gsl_vector_get(x,0);
    double mD = gsl_vector_get(x,1);
    double mS = gsl_vector_get(x,2);
    double effCP = gsl_vector_get(x,3);


    //define parameters
    SU3NJL3DCutoffEqualChemPotFixedTempRhoB solution(auxiliar);

    NJLDimensionfulCouplings couplings = solution.getParametersNJL().getDimensionfulCouplings();

    NJL3DCutoffRegularizationScheme reguScheme = solution.getParametersNJL().getNJL3DCutoffRegularizationScheme();
    double cutoff = solution.getParametersNJL().getThreeMomentumCutoff();

    double Nc = solution.getParametersNJL().getNumberOfColours();

    double sigmaIntegralPrecision = solution.getParametersNJL().getSigmaIntegralPrecision();
    double thermoIntegralPrecision = solution.getParametersNJL().getThermoIntegralPrecision();

    double m0U = solution.getParametersNJL().getUpQuarkCurrentMass();
    double m0D = solution.getParametersNJL().getDownQuarkCurrentMass();
    double m0S = solution.getParametersNJL().getStrangeQuarkCurrentMass();

    double T = solution.getTemperature();
    double rhoB = solution.getBaryonDensity();


    //Calculate the sigma field and density for each quark flavour
    double sigmaU = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCP, mU, sigmaIntegralPrecision);
    double sigmaD = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCP, mD, sigmaIntegralPrecision);
    double sigmaS = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCP, mS, sigmaIntegralPrecision);

    double rhoU = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCP, mU, thermoIntegralPrecision);
    double rhoD = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCP, mD, thermoIntegralPrecision);
    double rhoS = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCP, mS, thermoIntegralPrecision);

	//system of equations
    double f0 = SU3NJLNulledGapEquation(couplings, mU-m0U, sigmaU, sigmaD, sigmaS, rhoU, rhoD, rhoS);
    double f1 = SU3NJLNulledGapEquation(couplings, mD-m0D, sigmaD, sigmaS, sigmaU, rhoD, rhoS, rhoU);
    double f2 = SU3NJLNulledGapEquation(couplings, mS-m0S, sigmaS, sigmaU, sigmaD, rhoS, rhoU, rhoD);
	double f3 = rhoB - SU3BaryonDensity(rhoU, rhoD, rhoS);

	gsl_vector_set (f, 0, f0);
	gsl_vector_set (f, 1, f1);
	gsl_vector_set (f, 2, f2);
	gsl_vector_set (f, 3, f3);

	return GSL_SUCCESS;
}


bool SU3NJL3DCutoffEqualChemPotFixedTempRhoB::testSolution(double precision)
{   
    double x[4];
    x[0] = getUpQuarkEffectiveMass();
    x[1] = getDownQuarkEffectiveMass();
    x[2] = getStrangeQuarkEffectiveMass();
    x[3] = getQuarkEffectiveChemicalPotential();

    //the test below (gsl) resturns 0 if the sum_i abs(residual_i) < precision
    int gslTest = multiDimensionalRootFindTestResidual(4, precision, &x[0], this, &SU3NJL3DCutoffGapEquationsEqualChemicalPotentialFixedTemperatureBaryonDensity);

    if (gslTest==0){ return true; }
    else{ return false; }
}


double SU3NJL3DCutoffEqualChemPotFixedTempRhoB::calculatePressure(double vacuumPressure)
{
    double pressure = SU3NJL3DCutoffPressure(parametersNJL, temperature,
                                             upQuarkEffectiveMass, 
                                             downQuarkEffectiveMass, 
                                             strangeQuarkEffectiveMass, 
                                             quarkEffectiveChemicalPotential, 
                                             quarkEffectiveChemicalPotential, 
                                             quarkEffectiveChemicalPotential);

    pressure = pressure - vacuumPressure;

    return pressure;
}


double SU3NJL3DCutoffEqualChemPotFixedTempRhoB::calculateEnergyDensity(double vacuumEnergyDensity)
{
    double energy = SU3NJL3DCutoffEnergyDensity(parametersNJL, temperature,
                                                upQuarkEffectiveMass, 
                                                downQuarkEffectiveMass, 
                                                strangeQuarkEffectiveMass, 
                                                quarkEffectiveChemicalPotential, 
                                                quarkEffectiveChemicalPotential, 
                                                quarkEffectiveChemicalPotential);

    energy = energy - vacuumEnergyDensity;

    return energy;
}


double SU3NJL3DCutoffEqualChemPotFixedTempRhoB::calculateEntropyDensity()
{
    double entropy = SU3NJL3DCutoffEntropyDensity(parametersNJL, temperature,
                                                  upQuarkEffectiveMass, 
                                                  downQuarkEffectiveMass, 
                                                  strangeQuarkEffectiveMass, 
                                                  quarkEffectiveChemicalPotential, 
                                                  quarkEffectiveChemicalPotential, 
                                                  quarkEffectiveChemicalPotential);

    return entropy;
}


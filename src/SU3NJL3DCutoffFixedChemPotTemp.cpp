#include <cmath>
#include <iostream>
#include "rootSolverGSL.h"
#include "OneFermionLineIntegral.h"
#include "SU3NJL3DCutoff.h"
#include "SU3NJL3DCutoffFixedChemPotTemp.h"
#include "SU3NJL3DCutoffIntegratedCrossSections.h"

using namespace std;


SU3NJL3DCutoffFixedChemPotTemp::SU3NJL3DCutoffFixedChemPotTemp(void* auxiliar)
{	
	parametersNJL = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->parametersNJL;

	temperature = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->temperature;

	upQuarkChemicalPotential = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->upQuarkChemicalPotential;
	downQuarkChemicalPotential = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->downQuarkChemicalPotential;
	strangeQuarkChemicalPotential = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->strangeQuarkChemicalPotential;

	upQuarkEffectiveMass = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->upQuarkEffectiveMass;
	downQuarkEffectiveMass = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->downQuarkEffectiveMass;
	strangeQuarkEffectiveMass = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->strangeQuarkEffectiveMass;
};


SU3NJL3DCutoffFixedChemPotTemp::SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffParameters parametersNJLAux, 
															   double temperatureAux, 
															   double upQuarkChemicalPotentialAux,
															   double downQuarkChemicalPotentialAux,
															   double strangeQuarkChemicalPotentialAux)
{
	parametersNJL = parametersNJLAux;

	temperature = temperatureAux;

	upQuarkChemicalPotential = upQuarkChemicalPotentialAux;
	downQuarkChemicalPotential = downQuarkChemicalPotentialAux;
	strangeQuarkChemicalPotential = strangeQuarkChemicalPotentialAux;
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
    double cPU = solution.getUpQuarkChemicalPotential();
    double cPD = solution.getDownQuarkChemicalPotential();
    double cPS = solution.getStrangeQuarkChemicalPotential();

    //This solution does not take into account vector degrees of freedom
    lagrangianInteractions lagrangianInteractionsAux = solution.getParametersNJL().getDimensionfulCouplings().getLagrangianInteractions();
    if ( lagrangianInteractionsAux!=interactions_4SP_det && lagrangianInteractionsAux!=interactions_4SP_det_8SP  )
    {   
        cout << "Lagrangian interactions contain vector degrees of freedom! The class SU3NJL3DCutoffFixedChemPotTemp is not prepared for this!\n";
        abort();
    }

    //Calculate the sigma field and density for each quark flavour
    double sigmaU = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, cPU, mU, sigmaIntegralPrecision);
    double sigmaD = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, cPD, mD, sigmaIntegralPrecision);
    double sigmaS = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, cPS, mS, sigmaIntegralPrecision);

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



SU3NJL3DCutoffMeson SU3NJL3DCutoffFixedChemPotTemp::calculateMesonMassAndWidth(mesonState meson, double precision, MultiRootFindingMethod method, double mesonMassGuess, double mesonWidthGuess)
{   
    double mesonPropagatorPrecision = parametersNJL.getSigmaIntegralPrecision();

    SU3NJL3DCutoffMeson mesonAux(parametersNJL, temperature, 
                                 upQuarkChemicalPotential, downQuarkChemicalPotential, strangeQuarkChemicalPotential, 
                                 upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                                 mesonPropagatorPrecision, meson);

    mesonAux.calculateMesonMassAndWidth(precision, method, mesonMassGuess, mesonWidthGuess);

    return mesonAux;
}


std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(SU3NJL3DCutoffVacuum vacuumSol, double maxTemperature, int numberOfPoints, double precision, MultiRootFindingMethod method)
{
    double chemPotU = 0.0;
    double chemPotD = 0.0;
    double chemPotS = 0.0;

    double effMassU = vacuumSol.getUpQuarkEffectiveMass();
    double effMassD = vacuumSol.getDownQuarkEffectiveMass();
    double effMassS = vacuumSol.getStrangeQuarkEffectiveMass();

    double minTemperature = 1E-5;
    double deltaTemperature = (maxTemperature-minTemperature)/(numberOfPoints - 1);

    vector<SU3NJL3DCutoffFixedChemPotTemp> solutions;
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        double T = minTemperature + i*deltaTemperature;

        SU3NJL3DCutoffFixedChemPotTemp inMediumSol(vacuumSol.getParametersNJL(), T, chemPotU, chemPotD, chemPotS);
        inMediumSol.solve(precision, method, effMassU, effMassD, effMassS);

        effMassU = inMediumSol.getUpQuarkEffectiveMass();
        effMassD = inMediumSol.getDownQuarkEffectiveMass();
        effMassS = inMediumSol.getStrangeQuarkEffectiveMass();

        //cout << "Mu=" << effMassU << "GeV" << "\t" << "Md=" << effMassD << "GeV" << "\t" << "Ms=" << effMassS << "GeV" << "\n";

        if (inMediumSol.testSolution(1E-8)==true){ solutions.push_back(inMediumSol); }
    }

    return solutions;
}


std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromFiniteTemperatureToFiniteChemicalPotential(SU3NJL3DCutoffFixedChemPotTemp finiteTSol, double maxChemPot, int numberOfPoints, double precision, MultiRootFindingMethod method)
{
    double effMassU = finiteTSol.getUpQuarkEffectiveMass();
    double effMassD = finiteTSol.getDownQuarkEffectiveMass();
    double effMassS = finiteTSol.getStrangeQuarkEffectiveMass();

    double minChemPot = 1E-5;
    double delta = (maxChemPot-minChemPot)/(numberOfPoints - 1);

    vector<SU3NJL3DCutoffFixedChemPotTemp> solutions;
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        double chemPot = minChemPot + i*delta;

        SU3NJL3DCutoffFixedChemPotTemp inMediumSol(finiteTSol.getParametersNJL(), finiteTSol.getTemperature(), chemPot, chemPot, chemPot);
        inMediumSol.solve(precision, method, effMassU, effMassD, effMassS);

        effMassU = inMediumSol.getUpQuarkEffectiveMass();
        effMassD = inMediumSol.getDownQuarkEffectiveMass();
        effMassS = inMediumSol.getStrangeQuarkEffectiveMass();

        //cout << "Mu=" << effMassU << "GeV" << "\t" << "Md=" << effMassD << "GeV" << "\t" << "Ms=" << effMassS << "GeV" << "\n";

        if (inMediumSol.testSolution(1E-8)==true){ solutions.push_back(inMediumSol); }
    }

    return solutions;
}


vector<SU3NJL3DCutoffMeson> mesonPropertiesFromVacuumToFiniteTemperatureAtZeroChemicalPotential
(SU3NJL3DCutoffVacuum vacuumSolution, vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSolution, mesonState meson, double mesonPropertiesPrecision, MultiRootFindingMethod method, double mesonMassVacuumGuess, double mesonWidthVacuumGuess)
{
    //calculate meson mass and width in the vacuum and add it to vector
    SU3NJL3DCutoffMeson mesonVacuum = vacuumSolution.calculateMesonMassAndWidth(meson, mesonPropertiesPrecision, method, mesonMassVacuumGuess, mesonWidthVacuumGuess);

    vector<SU3NJL3DCutoffMeson> mesonFiniteT;
    mesonFiniteT.push_back( mesonVacuum );

    for (int i = 0; i < int(finiteTempSolution.size()); ++i)
    {   
        double mesonMassGuess, mesonWidthGuess;
        if ( i==0 )
        {   
            //use vacuum point as guess
            mesonMassGuess = mesonFiniteT[i].getMesonMass();
            mesonWidthGuess = mesonFiniteT[i].getMesonWidth();
        }
        else
        {   
            //use 2 previous points to find a guess
            double x1, y1, x2, y2, x;
            x1 = mesonFiniteT[i].getTemperature();
            x2 = mesonFiniteT[i-1].getTemperature();
            x = finiteTempSolution[i].getTemperature();

            y1 = mesonFiniteT[i].getMesonMass();
            y2 = mesonFiniteT[i-1].getMesonMass();
            mesonMassGuess = linearFit(x, x1, y1, x2, y2);

            y1 = mesonFiniteT[i].getMesonWidth();
            y2 = mesonFiniteT[i-1].getMesonWidth();
            mesonWidthGuess = linearFit(x, x1, y1, x2, y2);
        }
        
        SU3NJL3DCutoffMeson mesonFiniteTAux = finiteTempSolution[i].calculateMesonMassAndWidth(meson, mesonPropertiesPrecision, method, mesonMassGuess, mesonWidthGuess);
        mesonFiniteT.push_back( mesonFiniteTAux );
    }

    return mesonFiniteT;
}



//Evaluate Cross section for paper at finite chemical potential using Klevansky parameter set
void evaluateCrossSectionsPaperWithKlevanskyParameterSet(double T, double chemPot, int numberOfCrossSectionPoints)
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


    //numerical precisions
    double gapPrecision = 1E-8;
    double mesonPropagatorIntegralPrecision = 1E-8;
    double crossSectionIntegralPrecision = 1E-4;


    //find solution in the vacuum for the Klevansky parameter set
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, hybrids, 0.3, 0.3, 0.5);

    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";


    //solve gap equation from the vacuum up to finite temperature
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSol = solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, T, 100, gapPrecision, hybrids);

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
        vector<SU3NJL3DCutoffFixedChemPotTemp> inMediumSol = solveFromFiniteTemperatureToFiniteChemicalPotential(finiteTempSol[int(finiteTempSol.size()-1)], chemPot, 100, gapPrecision, hybrids);

        effMassU = inMediumSol[int(inMediumSol.size()-1)].getUpQuarkEffectiveMass();
        effMassD = inMediumSol[int(inMediumSol.size()-1)].getDownQuarkEffectiveMass();
        effMassS = inMediumSol[int(inMediumSol.size()-1)].getStrangeQuarkEffectiveMass();

        cout << "Mu=" << effMassU << "GeV" << "\t" 
             << "Md=" << effMassD << "GeV" << "\t" 
             << "Ms=" << effMassS << "GeV" << "\n";
    }


    //evaluate cross sections
    evaluateCrossSectionsPaperFiniteChemicalPotential(parameters, T, 
                                                      chemPot, chemPot, chemPot, 
                                                      effMassU, effMassD, effMassS, 
                                                      mesonPropagatorIntegralPrecision,
                                                      false, crossSectionIntegralPrecision,
                                                      numberOfCrossSectionPoints);
}



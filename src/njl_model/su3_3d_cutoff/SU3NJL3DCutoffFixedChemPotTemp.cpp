#include <cmath>
#include <iostream>
#include "gsl_wrapper/root_solver_gsl.h"
#include "njl_model/n_fermion_line_integrals/one_fermion_line_integral_3d_cutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.h"

using namespace std;


//linear interpolation and extrapolation
double linearFit(double x, double x1, double y1, double x2, double y2)
{
    double m = ( y1 - y2 )/( x1 - x2 );
    double b = 0.5*( ( y1 + y2 ) - m*( x1 + x2 ) );
    double y = m*x + b;

    return y;
}


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
    mesonID = ((class SU3NJL3DCutoffFixedChemPotTemp *)(auxiliar))->mesonID;
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


SU3NJL3DCutoffFixedChemPotTemp::SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffParameters parametersNJLAux, 
                                                               double upQuarkChemicalPotentialAux,
                                                               double downQuarkChemicalPotentialAux,
                                                               double strangeQuarkChemicalPotentialAux)
{
    parametersNJL = parametersNJLAux;

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
    LagrangianInteractions lagrangianInteractionsAux = solution.getParametersNJL().getDimensionfulCouplings().getLagrangianInteractions();
    if ( lagrangianInteractionsAux!=SP4Q_DET2NFQ && lagrangianInteractionsAux!=SP4Q_DET2NFQ_SP8Q  )
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



SU3NJL3DCutoffMeson SU3NJL3DCutoffFixedChemPotTemp::calculateMesonMassAndWidth(mesonState mesonIDAux, double precision, MultiRootFindingMethod method, double mesonMassGuess, double mesonWidthGuess)
{   
    double mesonPropagatorPrecision = parametersNJL.getSigmaIntegralPrecision();

    SU3NJL3DCutoffMeson mesonAux(parametersNJL, temperature, 
                                 upQuarkChemicalPotential, downQuarkChemicalPotential, strangeQuarkChemicalPotential, 
                                 upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                                 mesonPropagatorPrecision, mesonIDAux);

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
        if (inMediumSol.testSolution(precision)==true){ solutions.push_back(inMediumSol); }
    }

    return solutions;
}


std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromLowToHighTemperatureAtZeroChemicalPotential(SU3NJL3DCutoffFixedChemPotTemp minTemperatureSolution, double maxTemperature, int numberOfPoints, double precision, MultiRootFindingMethod method)
{
    double chemPotU = 0.0;
    double chemPotD = 0.0;
    double chemPotS = 0.0;

    double effMassU = minTemperatureSolution.getUpQuarkEffectiveMass();
    double effMassD = minTemperatureSolution.getDownQuarkEffectiveMass();
    double effMassS = minTemperatureSolution.getStrangeQuarkEffectiveMass();

    double minTemperature = minTemperatureSolution.getTemperature();
    double deltaTemperature = (maxTemperature-minTemperature)/(numberOfPoints - 1);

    vector<SU3NJL3DCutoffFixedChemPotTemp> solutions;
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        double T = minTemperature + i*deltaTemperature;

        SU3NJL3DCutoffFixedChemPotTemp inMediumSol(minTemperatureSolution.getParametersNJL(), T, chemPotU, chemPotD, chemPotS);
        inMediumSol.solve(precision, method, effMassU, effMassD, effMassS);

        effMassU = inMediumSol.getUpQuarkEffectiveMass();
        effMassD = inMediumSol.getDownQuarkEffectiveMass();
        effMassS = inMediumSol.getStrangeQuarkEffectiveMass();

        //cout << "Mu=" << effMassU << "GeV" << "\t" << "Md=" << effMassD << "GeV" << "\t" << "Ms=" << effMassS << "GeV" << "\n";
        if (inMediumSol.testSolution(precision)==true){ solutions.push_back(inMediumSol); }
    }

    return solutions;
}


std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromLowToHighTemperature(SU3NJL3DCutoffFixedChemPotTemp minTemperatureSolution, double maxTemperature, int numberOfPoints, double precision, MultiRootFindingMethod method)
{
    double chemPotU = minTemperatureSolution.getUpQuarkChemicalPotential();
    double chemPotD = minTemperatureSolution.getDownQuarkChemicalPotential();
    double chemPotS = minTemperatureSolution.getStrangeQuarkChemicalPotential();

    double effMassU = minTemperatureSolution.getUpQuarkEffectiveMass();
    double effMassD = minTemperatureSolution.getDownQuarkEffectiveMass();
    double effMassS = minTemperatureSolution.getStrangeQuarkEffectiveMass();

    double minTemperature = minTemperatureSolution.getTemperature();
    double deltaTemperature = (maxTemperature-minTemperature)/(numberOfPoints - 1);

    vector<SU3NJL3DCutoffFixedChemPotTemp> solutions;
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        double T = minTemperature + i*deltaTemperature;

        SU3NJL3DCutoffFixedChemPotTemp inMediumSol(minTemperatureSolution.getParametersNJL(), T, chemPotU, chemPotD, chemPotS);
        inMediumSol.solve(precision, method, effMassU, effMassD, effMassS);

        effMassU = inMediumSol.getUpQuarkEffectiveMass();
        effMassD = inMediumSol.getDownQuarkEffectiveMass();
        effMassS = inMediumSol.getStrangeQuarkEffectiveMass();

        //cout << "Mu=" << effMassU << "GeV" << "\t" << "Md=" << effMassD << "GeV" << "\t" << "Ms=" << effMassS << "GeV" << "\n";
        if (inMediumSol.testSolution(precision)==true){ solutions.push_back(inMediumSol); }
    }

    return solutions;
}


vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromFiniteTemperatureToFiniteChemicalPotential(SU3NJL3DCutoffFixedChemPotTemp finiteTSol, double maxChemPot, int numberOfPoints, double precision, MultiRootFindingMethod method)
{
    double effMassU = finiteTSol.getUpQuarkEffectiveMass();
    double effMassD = finiteTSol.getDownQuarkEffectiveMass();
    double effMassS = finiteTSol.getStrangeQuarkEffectiveMass();

    double minChemPot = finiteTSol.getUpQuarkChemicalPotential();
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

        if (inMediumSol.testSolution(precision)==true){ solutions.push_back(inMediumSol); }
    }

    return solutions;
}


vector<SU3NJL3DCutoffMeson> mesonPropertiesFromVacuumToFiniteTemperatureAtZeroChemicalPotential
(SU3NJL3DCutoffVacuum vacuumSolution, vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSolution, mesonState mesonID, double mesonPropertiesPrecision, MultiRootFindingMethod method, double mesonMassVacuumGuess, double mesonWidthVacuumGuess)
{
    //calculate meson mass and width in the vacuum and add it to vector
    SU3NJL3DCutoffMeson mesonVacuum = vacuumSolution.calculateMesonMassAndWidth(mesonID, mesonPropertiesPrecision, method, mesonMassVacuumGuess, mesonWidthVacuumGuess);

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
        
        SU3NJL3DCutoffMeson mesonFiniteTAux = finiteTempSolution[i].calculateMesonMassAndWidth(mesonID, mesonPropertiesPrecision, method, mesonMassGuess, mesonWidthGuess);
        mesonFiniteT.push_back( mesonFiniteTAux );
    }

    return mesonFiniteT;
}


int SU3NJL3DCutoffNondiagonalMesonMottTemperatureFixedChemicalPotentials(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{
    //define variables
    double mU = gsl_vector_get(x,0);
    double mD = gsl_vector_get(x,1);
    double mS = gsl_vector_get(x,2);
    double mottTemp = gsl_vector_get(x,3);

    //Gap equations stuff
    SU3NJL3DCutoffFixedChemPotTemp solution(auxiliar);
    solution.setTemperature( mottTemp );
    
    //Meson propagator stuff
    double mesonPropagatorPrecision = solution.getParametersNJL().getSigmaIntegralPrecision();
    mesonState mesonID = solution.getMesonID();
    double k0 = mesonStateMassAtMeltingPoint(mU, mD, mS, mesonID);
    gsl_complex mesonPropagator =
    nonDiagonalMesonPropagator(solution.getParametersNJL(), 
                               mottTemp, 
                               solution.getUpQuarkChemicalPotential(), 
                               solution.getDownQuarkChemicalPotential(), 
                               solution.getStrangeQuarkChemicalPotential(), 
                               mU, mD, mS, k0, 0.0, 0.0, mesonPropagatorPrecision, mesonID);

    //Systems of equations to solve
    SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature(x, &solution, &f[0]); //use gap equations for this physical scenario
    double f0 = gsl_vector_get(f,0);
    double f1 = gsl_vector_get(f,1);
    double f2 = gsl_vector_get(f,2);
    double f3 = GSL_REAL( gsl_complex_inverse( mesonPropagator ) );

    gsl_vector_set (f, 0, f0);
    gsl_vector_set (f, 1, f1);
    gsl_vector_set (f, 2, f2);
    gsl_vector_set (f, 3, f3);

    return GSL_SUCCESS;
}


void SU3NJL3DCutoffFixedChemPotTemp::findNondiagonalMesonMottTemperature(mesonState mesonIDAux,
                                                                         double precision, 
                                                                         MultiRootFindingMethod method, 
                                                                         double upQuarkEffectiveMassGuess, 
                                                                         double downQuarkEffectiveMassGuess, 
                                                                         double strangeQuarkEffectiveMassGuess,
                                                                         double mottTemperatureGuess)
{   
    double x[4];
    x[0] = upQuarkEffectiveMassGuess; 
    x[1] = downQuarkEffectiveMassGuess; 
    x[2] = strangeQuarkEffectiveMassGuess;
    x[3] = mottTemperatureGuess;
    
    setMesonID(mesonIDAux);
    multiDimensionalRootFind(4, precision, &x[0], this, &SU3NJL3DCutoffNondiagonalMesonMottTemperatureFixedChemicalPotentials, method);

    setUpQuarkEffectiveMass(x[0]);
    setDownQuarkEffectiveMass(x[1]);
    setStrangeQuarkEffectiveMass(x[2]);
    setTemperature(x[3]);
}


SU3NJL3DCutoffFixedChemPotTemp nondiagonalMesonMeltingPoint(SU3NJL3DCutoffVacuum vacuumSolution, vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSolution, mesonState mesonID, double mesonPropertiesPrecision, MultiRootFindingMethod method, double mesonMassVacuumGuess, double mesonWidthVacuumGuess)
{   
    //calculate the meson mass behaviour for the provided range of NJL solutions
    vector<SU3NJL3DCutoffMeson> mesonFiniteT = 
    mesonPropertiesFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuumSolution, finiteTempSolution, mesonID, mesonPropertiesPrecision, method, mesonMassVacuumGuess, mesonWidthVacuumGuess);


    //calculate Mott temperature of a given meson if it exists
    int meltingPointGuess = 0;
    double meltTest = 0;
    for (int i = 0; i < int(mesonFiniteT.size()); ++i)
    {   
        double k0AtMeltingPoint = mesonStateMassAtMeltingPoint(mesonFiniteT[i].getUpQuarkEffectiveMass(), 
                                                               mesonFiniteT[i].getDownQuarkEffectiveMass(), 
                                                               mesonFiniteT[i].getStrangeQuarkEffectiveMass(), 
                                                               mesonID);
        meltTest = mesonFiniteT[i].getMesonMass() - k0AtMeltingPoint;
        if ( meltTest>0 )
        { 
            meltingPointGuess = i; 
            break; 
        }
    }

    SU3NJL3DCutoffFixedChemPotTemp mottSolution;
    //If there is apoint for which the meson mass is approximately equal the sum of two quark masses, there is a melting point
    if ( meltTest>0 )
    {
        mottSolution = SU3NJL3DCutoffFixedChemPotTemp(vacuumSolution.getParametersNJL(), 0.0, 0.0, 0.0);
        mottSolution.findNondiagonalMesonMottTemperature(mesonID, mesonPropertiesPrecision, HYBRIDS, 
                                                         mesonFiniteT[meltingPointGuess].getUpQuarkEffectiveMass(), 
                                                         mesonFiniteT[meltingPointGuess].getDownQuarkEffectiveMass(), 
                                                         mesonFiniteT[meltingPointGuess].getStrangeQuarkEffectiveMass(), 
                                                         mesonFiniteT[meltingPointGuess].getTemperature());
        cout << mottSolution.getTemperature() << "\n";
    }

    return mottSolution;
}


//Evaluate Cross section for paper at finite chemical potential using Klevansky parameter set
void evaluateCrossSectionsPaperWithKlevanskyParameterSet(double T, double chemPot, int numberOfCrossSectionPoints, int numberOfThreads)
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
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);


    //numerical precisions
    double gapPrecision = 1E-8;
    double mesonPropagatorIntegralPrecision = 1E-8;
    double crossSectionIntegralPrecision = 1E-4;


    //find solution in the vacuum for the Klevansky parameter set
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, HYBRIDS, 0.3, 0.3, 0.5);

    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";


    //solve gap equation from the vacuum up to finite temperature
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSol = solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, T, 100, gapPrecision, HYBRIDS);

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
        vector<SU3NJL3DCutoffFixedChemPotTemp> inMediumSol = solveFromFiniteTemperatureToFiniteChemicalPotential(finiteTempSol[int(finiteTempSol.size()-1)], chemPot, 100, gapPrecision, HYBRIDS);

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
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ, gs, kappa);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setA");


    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, HYBRIDS, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";


    cout << "Pseudoscalar meson masses:\n";

    SU3NJL3DCutoffMeson pionPlusMassSolution = vacuum.calculateMesonMassAndWidth(pionPlus, 1E-7, HYBRIDS, 0.2, 0.2);
    cout << pionPlusMassSolution.getMesonMass() << "\t" << pionPlusMassSolution.getMesonWidth() << "\n";

    SU3NJL3DCutoffMeson kaonPlusMassSolution = vacuum.calculateMesonMassAndWidth(kaonPlus, 1E-7, HYBRIDS, 0.4, 0.2);
    cout << kaonPlusMassSolution.getMesonMass() << "\t" << kaonPlusMassSolution.getMesonWidth() << "\n";

    SU3NJL3DCutoffMeson diagonalMassSolution;
    diagonalMassSolution = vacuum.calculateMesonMassAndWidth(diagonalPseudoscalars, 1E-7, HYBRIDS, 0.5, 0.2);
    cout << diagonalMassSolution.getMesonMass() << "\t" << diagonalMassSolution.getMesonWidth() << "\n";

    diagonalMassSolution = vacuum.calculateMesonMassAndWidth(diagonalPseudoscalars, 1E-7, HYBRIDS, 1.0, 0.2);
    cout << diagonalMassSolution.getMesonMass() << "\t" << diagonalMassSolution.getMesonWidth() << "\n";

     cout << "Scalar meson masses:\n";

    SU3NJL3DCutoffMeson sigmaPionPlusMassSolution = vacuum.calculateMesonMassAndWidth(sigmaPionPlus, 1E-7, HYBRIDS, 0.2, 0.2);
    cout << sigmaPionPlusMassSolution.getMesonMass() << "\t" << sigmaPionPlusMassSolution.getMesonWidth() << "\n";


    //solve model at zero chemical potential up to some finite temperature
    double maximumTemperature = 0.400;
    int numberOfPoints = 400;
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTSolution = 
    solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, maximumTemperature, numberOfPoints, gapPrecision, HYBRIDS);


    //Search for meson melting points in the range of temperatures considered above
    double mesonPropertiesPrecision = 1E-7;
    double mesonMassVacuumGuess;
    double mesonWidthVacuumGuess;
    mesonState mesonID;

    mesonMassVacuumGuess = 0.2;
    mesonWidthVacuumGuess = 0.2;
    mesonID = pionPlus;
    SU3NJL3DCutoffFixedChemPotTemp meltingPointPionPlus = nondiagonalMesonMeltingPoint(vacuum, finiteTSolution, mesonID, mesonPropertiesPrecision, HYBRIDS, mesonMassVacuumGuess, mesonWidthVacuumGuess);
    cout << "pionPlus TMott: " << meltingPointPionPlus.getTemperature() << "\n";

    mesonMassVacuumGuess = 0.5;
    mesonWidthVacuumGuess = 0.2;
    mesonID = kaonPlus;
    SU3NJL3DCutoffFixedChemPotTemp meltingPointKaonPlus = nondiagonalMesonMeltingPoint(vacuum, finiteTSolution, mesonID, mesonPropertiesPrecision, HYBRIDS, mesonMassVacuumGuess, mesonWidthVacuumGuess);
    cout << "kaonPlus TMott: " << meltingPointKaonPlus.getTemperature() << "\n";
}


void SU3NJL3DCutoffFixedChemPotTemp::evaluateCrossSectionsEqualLightMasses(
    SU3NJL3DCutoffParameters& parameters,                                    
    double precisionVacuum,                                    
    MultiRootFindingMethod methodVacuum,                                    
    double lightQuarkMassGuess, 
    double strangeQuarkMassGuess,
    double temperature,
    int numberOfPointsFromVacToFinTemp,
    double precisionVacToFinTemp,
    MultiRootFindingMethod methodVacToFinTemp,
    double quarkChemicalPotential,
    int numberOfPointsFromFinTempToFinChemPot,
    double precisionFinTempToFinChemPot,
    MultiRootFindingMethod methodFinTempToFinChemPot,
    double propagatorIntegralPrecision,
    bool largeAngleScatteringContribution,
    double precisionCrossSections,
    int numberOfPointsCrossSections,
    int numberOfThreads
)
{
    SU3NJL3DCutoffVacuum vacuum = SU3NJL3DCutoffVacuum::calculateVacuumMasses(
        parameters,                                    
        precisionVacuum,                                    
        methodVacuum,                                    
        lightQuarkMassGuess, 
        lightQuarkMassGuess, 
        strangeQuarkMassGuess
    );

    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSol = solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
        vacuum, 
        temperature, 
        numberOfPointsFromVacToFinTemp, 
        precisionVacToFinTemp, 
        methodVacToFinTemp
    );
    
    double effMassU = finiteTempSol[int(finiteTempSol.size()-1)].getUpQuarkEffectiveMass();
    double effMassD = finiteTempSol[int(finiteTempSol.size()-1)].getDownQuarkEffectiveMass();
    double effMassS = finiteTempSol[int(finiteTempSol.size()-1)].getStrangeQuarkEffectiveMass();

    cout << "\n";
    cout << "Temperature[GeV] = " << temperature << "\n";
    cout << "quarkChemicalPotential[GeV] = 0.0 \n";
    cout << "Effective masses: \n";
    cout << "Mu[GeV] = " << effMassU << "\n" 
         << "Md[GeV] = " << effMassD << "\n" 
         << "Ms[GeV] = " << effMassS << "\n";

    if ( quarkChemicalPotential>0.0 )
    {
        vector<SU3NJL3DCutoffFixedChemPotTemp> inMediumSol = solveFromFiniteTemperatureToFiniteChemicalPotential(
            finiteTempSol[int(finiteTempSol.size()-1)], 
            quarkChemicalPotential, 
            numberOfPointsFromFinTempToFinChemPot, 
            precisionFinTempToFinChemPot, 
            methodFinTempToFinChemPot
        );

        effMassU = inMediumSol[int(inMediumSol.size()-1)].getUpQuarkEffectiveMass();
        effMassD = inMediumSol[int(inMediumSol.size()-1)].getDownQuarkEffectiveMass();
        effMassS = inMediumSol[int(inMediumSol.size()-1)].getStrangeQuarkEffectiveMass();

        cout << "\n";
        cout << "Temperature[GeV] = " << temperature << "\n";
        cout << "quarkChemicalPotential[GeV] = " << quarkChemicalPotential << "\n";
        cout << "Effective masses: \n";
        cout << "Mu[GeV] = " << effMassU << "\n" 
             << "Md[GeV] = " << effMassD << "\n" 
             << "Ms[GeV] = " << effMassS << "\n";
    }

    evaluateCrossSectionsEqualLightMassesEqualChemicalPotential(
        parameters, temperature, 
        quarkChemicalPotential, quarkChemicalPotential, quarkChemicalPotential, 
        effMassU, effMassD, effMassS, 
        propagatorIntegralPrecision,
        largeAngleScatteringContribution, 
        precisionCrossSections,
        numberOfPointsCrossSections,
        numberOfThreads
    );
}

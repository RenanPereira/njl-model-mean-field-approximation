#include <cmath>
#include <iostream>
#include "gsl_wrapper/root_solver_gsl.h"
#include "njl_model/n_fermion_line_integrals/one_fermion_line_integral_3d_cutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.h"
#include "math_utils/useful_functions.h"
#include "utils/format_utils.h"

using namespace std;

// Constructor that takes generic pointer (void*), casts it into class object and copy it to the current instance
SU3NJL3DCutoffFixedChemPotTemp::SU3NJL3DCutoffFixedChemPotTemp(void* auxiliar)
{
    SU3NJL3DCutoffFixedChemPotTemp* solution = static_cast<SU3NJL3DCutoffFixedChemPotTemp*>(auxiliar);
    *this = *solution;
}

SU3NJL3DCutoffFixedChemPotTemp::SU3NJL3DCutoffFixedChemPotTemp(
    SU3NJL3DCutoffParameters parametersNJLAux, 
    double temperatureAux, 
    double upQuarkChemicalPotentialAux, 
    double downQuarkChemicalPotentialAux, 
    double strangeQuarkChemicalPotentialAux
)
{
	parametersNJL = parametersNJLAux;

	temperature = temperatureAux;
	upQuarkChemicalPotential = upQuarkChemicalPotentialAux;
	downQuarkChemicalPotential = downQuarkChemicalPotentialAux;
	strangeQuarkChemicalPotential = strangeQuarkChemicalPotentialAux;
}

SU3NJL3DCutoffFixedChemPotTemp::SU3NJL3DCutoffFixedChemPotTemp(
    SU3NJL3DCutoffParameters parametersNJLAux, 
    double upQuarkChemicalPotentialAux, 
    double downQuarkChemicalPotentialAux, 
    double strangeQuarkChemicalPotentialAux 
)
{
    parametersNJL = parametersNJLAux;

    upQuarkChemicalPotential = upQuarkChemicalPotentialAux;
    downQuarkChemicalPotential = downQuarkChemicalPotentialAux;
    strangeQuarkChemicalPotential = strangeQuarkChemicalPotentialAux;
}

SU3NJL3DCutoffFixedChemPotTemp::SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffVacuum &vacuumSolution)
{
    parametersNJL = vacuumSolution.getParametersNJL();

    temperature = 0.0;
    upQuarkChemicalPotential = 0.0;
    downQuarkChemicalPotential = 0.0;
    strangeQuarkChemicalPotential = 0.0;
    
    upQuarkEffectiveMass = vacuumSolution.getUpQuarkEffectiveMass();
	downQuarkEffectiveMass = vacuumSolution.getDownQuarkEffectiveMass();
	strangeQuarkEffectiveMass = vacuumSolution.getStrangeQuarkEffectiveMass();

    pressure = 0.0;
	energyDensity = 0.0;
	entropyDensity = 0.0;
}

void SU3NJL3DCutoffFixedChemPotTemp::solve(
    double precision, 
    MultiRootFindingMethod method, 
    double upQuarkEffectiveMassGuess, 
    double downQuarkEffectiveMassGuess, 
    double strangeQuarkEffectiveMassGuess
)
{	
	double x[3];
    x[0] = upQuarkEffectiveMassGuess; 
    x[1] = downQuarkEffectiveMassGuess; 
    x[2] = strangeQuarkEffectiveMassGuess;
    
    multiDimensionalRootFind(
        3, 
        precision, 
        &x[0], 
        this, 
        &SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature, 
        method
    );

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
    
    //cast void* auxiliar into SU3NJL3DCutoffFixedChemPotTemp* and get necessary parameters
    const SU3NJL3DCutoffFixedChemPotTemp* solution = static_cast<const SU3NJL3DCutoffFixedChemPotTemp*>(auxiliar);

    NJLDimensionfulCouplings couplings = solution->getParametersNJL().getDimensionfulCouplings();
    LagrangianInteractions lagrangianInteractions = solution->getParametersNJL().getDimensionfulCouplings().getLagrangianInteractions();

    NJL3DCutoffRegularizationScheme reguScheme = solution->getParametersNJL().getNJL3DCutoffRegularizationScheme();
    double cutoff = solution->getParametersNJL().getThreeMomentumCutoff();

    double Nc = solution->getParametersNJL().getNumberOfColours();

    double sigmaIntegralPrecision = solution->getParametersNJL().getSigmaIntegralPrecision();

    double m0U = solution->getParametersNJL().getUpQuarkCurrentMass();
    double m0D = solution->getParametersNJL().getDownQuarkCurrentMass();
    double m0S = solution->getParametersNJL().getStrangeQuarkCurrentMass();

    double T = solution->getTemperature();
    double cPU = solution->getUpQuarkChemicalPotential();
    double cPD = solution->getDownQuarkChemicalPotential();
    double cPS = solution->getStrangeQuarkChemicalPotential();

    //This solution does not take into account vector degrees of freedom
    if ( lagrangianInteractions!=SP4Q_DET2NFQ && lagrangianInteractions!=SP4Q_DET2NFQ_SP8Q  )
    {   
        cout << "Lagrangian interactions contain vector degrees of freedom! "
             << "The class SU3NJL3DCutoffFixedChemPotTemp is not prepared for this!\n";
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
    int gslTest = multiDimensionalRootFindTestResidual(
        3, 
        precision, 
        &x[0], 
        this, 
        &SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature
    );

    if ( gslTest == 0 )
    { 
        return true; 
    }
    else
    { 
        return false; 
    }
}

SU3NJL3DCutoffMeson SU3NJL3DCutoffFixedChemPotTemp::calculateMesonMassAndWidth(
    mesonState mesonIDAux, 
    double precision, 
    MultiRootFindingMethod method, 
    double mesonMassGuess, 
    double mesonWidthGuess
)
{   
    double mesonPropagatorPrecision = parametersNJL.getSigmaIntegralPrecision();

    SU3NJL3DCutoffMeson mesonAux(
        parametersNJL, 
        temperature, 
        upQuarkChemicalPotential, 
        downQuarkChemicalPotential, 
        strangeQuarkChemicalPotential, 
        upQuarkEffectiveMass, 
        downQuarkEffectiveMass, 
        strangeQuarkEffectiveMass, 
        mesonPropagatorPrecision, 
        mesonIDAux
    );

    mesonAux.calculateMesonMassAndWidth(precision, method, mesonMassGuess, mesonWidthGuess);

    return mesonAux;
}

vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
    SU3NJL3DCutoffVacuum vacuumSol, 
    double nearVacuumTemperature,
    double maxTemperature, 
    int numberOfPoints, 
    double precision, 
    MultiRootFindingMethod method
)
{
    //Solve the model from the vacuum to some finite temperature at zero chemical potential
    cout << "\nSolving the SU3 NJL model, regularized by a 3D Cutoff, from the vacuum to some "
         << "finite temperature at zero chemical potential...\n";

    double chemPotU = 0.0;
    double chemPotD = 0.0;
    double chemPotS = 0.0;

    double effMassU = vacuumSol.getUpQuarkEffectiveMass();
    double effMassD = vacuumSol.getDownQuarkEffectiveMass();
    double effMassS = vacuumSol.getStrangeQuarkEffectiveMass();

    double deltaTemperature = (maxTemperature-nearVacuumTemperature)/(numberOfPoints - 1);

    vector<SU3NJL3DCutoffFixedChemPotTemp> solutions;
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        double T = nearVacuumTemperature + i*deltaTemperature;

        SU3NJL3DCutoffFixedChemPotTemp inMediumSol(vacuumSol.getParametersNJL(), T, chemPotU, chemPotD, chemPotS);
        inMediumSol.solve(precision, method, effMassU, effMassD, effMassS);

        effMassU = inMediumSol.getUpQuarkEffectiveMass();
        effMassD = inMediumSol.getDownQuarkEffectiveMass();
        effMassS = inMediumSol.getStrangeQuarkEffectiveMass();

        //cout << "Mu=" << effMassU << "GeV" << "\t" << "Md=" << effMassD << "GeV" << "\t" << "Ms=" << effMassS << "GeV" << "\n";
        
        //only store solutions which pass the test
        if (inMediumSol.testSolution(precision)==true)
        { 
            solutions.push_back(inMediumSol); 
        }
    }

    return solutions;
}

vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromLowToHighTemperatureAtZeroChemicalPotential(
    SU3NJL3DCutoffFixedChemPotTemp minTemperatureSolution, 
    double maxTemperature, 
    int numberOfPoints, 
    double precision, 
    MultiRootFindingMethod method
)
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
        if (inMediumSol.testSolution(precision)==true)
        { 
            solutions.push_back(inMediumSol); 
        }
    }

    return solutions;
}

vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromLowToHighTemperature(
    SU3NJL3DCutoffFixedChemPotTemp minTemperatureSolution, 
    double maxTemperature, 
    int numberOfPoints, 
    double precision, 
    MultiRootFindingMethod method
)
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
        if (inMediumSol.testSolution(precision)==true)
        { 
            solutions.push_back(inMediumSol); 
        }
    }

    return solutions;
}

vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromFiniteTemperatureToFiniteChemicalPotential(
    SU3NJL3DCutoffFixedChemPotTemp finiteTSol, 
    double maxChemPot, 
    int numberOfPoints, 
    double precision, 
    MultiRootFindingMethod method
)
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

vector<SU3NJL3DCutoffFixedChemPotTemp> solveVacuumToFiniteChemicalPotential(
    SU3NJL3DCutoffVacuum& vacuumSolution, 
    double maxChemPot, 
    int numberOfPoints, 
    double precision, 
    MultiRootFindingMethod method
)
{
    if (numberOfPoints < 2)
    {
        cout << "The parameter numberOfPoints must be >= 2 when calling function: " << __func__  << endl;
        abort();
    }

    SU3NJL3DCutoffParameters parameters = vacuumSolution.getParametersNJL();
    double effMassU = vacuumSolution.getUpQuarkEffectiveMass();
    double effMassD = vacuumSolution.getDownQuarkEffectiveMass();
    double effMassS = vacuumSolution.getStrangeQuarkEffectiveMass();

    double minChemPot = 0.0;
    double delta = (maxChemPot-minChemPot)/(numberOfPoints - 1);

    vector<SU3NJL3DCutoffFixedChemPotTemp> solutions;
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        double chemPot = minChemPot + i*delta;
        
        // create instance with temperature set to zero and equal chemical potentials
        SU3NJL3DCutoffFixedChemPotTemp inMediumSol(parameters, 0.0, chemPot, chemPot, chemPot);
        inMediumSol.solve(precision, method, effMassU, effMassD, effMassS);

        effMassU = inMediumSol.getUpQuarkEffectiveMass();
        effMassD = inMediumSol.getDownQuarkEffectiveMass();
        effMassS = inMediumSol.getStrangeQuarkEffectiveMass();

        //cout << "Mu=" << effMassU << "GeV" << "\t" << "Md=" << effMassD << "GeV" << "\t" << "Ms=" << effMassS << "GeV" << "\n";

        if (inMediumSol.testSolution(precision)==true)
        { 
            solutions.push_back(inMediumSol); 
        }
        else
        {
            cout << "Failed to reach desired precision when solving the model in function: "  << __func__  << endl;
            abort();
        }
    }

    return solutions;
}

vector<SU3NJL3DCutoffFixedChemPotTemp> solveUpToTemperature(
    SU3NJL3DCutoffFixedChemPotTemp solution, 
    double temperature, 
    int numberOfPoints, 
    double precision, 
    MultiRootFindingMethod method
)
{
    if (numberOfPoints < 2)
    {
        cout << "The parameter numberOfPoints must be >= 2 when calling function: " << __func__  << endl;
        abort();
    }

    SU3NJL3DCutoffParameters parameters = solution.getParametersNJL();
    double chemPotU = solution.getUpQuarkChemicalPotential();
    double chemPotD = solution.getDownQuarkChemicalPotential();
    double chemPotS = solution.getStrangeQuarkChemicalPotential();
    double effMassU = solution.getUpQuarkEffectiveMass();
    double effMassD = solution.getDownQuarkEffectiveMass();
    double effMassS = solution.getStrangeQuarkEffectiveMass();

    double minTemperature = solution.getTemperature();
    double deltaTemperature = (temperature-minTemperature)/(numberOfPoints - 1);

    vector<SU3NJL3DCutoffFixedChemPotTemp> solutions;
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        double T = minTemperature + i*deltaTemperature;

        SU3NJL3DCutoffFixedChemPotTemp inMediumSol(parameters, T, chemPotU, chemPotD, chemPotS);
        inMediumSol.solve(precision, method, effMassU, effMassD, effMassS);

        effMassU = inMediumSol.getUpQuarkEffectiveMass();
        effMassD = inMediumSol.getDownQuarkEffectiveMass();
        effMassS = inMediumSol.getStrangeQuarkEffectiveMass();

        //cout << "Mu=" << effMassU << "GeV" << "\t" << "Md=" << effMassD << "GeV" << "\t" << "Ms=" << effMassS << "GeV" << "\n";
        if (inMediumSol.testSolution(precision)==true)
        { 
            solutions.push_back(inMediumSol); 
        }
        else
        {
            cout << "Failed to reach desired precision when solving the model in function: "  << __func__  << endl;
            abort();
        }
    }

    return solutions;
}

vector<SU3NJL3DCutoffMeson> mesonPropertiesFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
    SU3NJL3DCutoffVacuum vacuumSolution, 
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSolution, 
    mesonState mesonID, 
    double mesonPropertiesPrecision, 
    MultiRootFindingMethod method, 
    double mesonMassVacuumGuess, 
    double mesonWidthVacuumGuess
)
{
    //calculate meson mass and width in the vacuum and add it to vector
    SU3NJL3DCutoffMeson mesonVacuum = vacuumSolution.calculateMesonMassAndWidth(
        mesonID, 
        mesonPropertiesPrecision, 
        method, 
        mesonMassVacuumGuess, 
        mesonWidthVacuumGuess
    );

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
        
        SU3NJL3DCutoffMeson mesonFiniteTAux = finiteTempSolution[i].calculateMesonMassAndWidth(
            mesonID, 
            mesonPropertiesPrecision, 
            method, 
            mesonMassGuess, 
            mesonWidthGuess
        );
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
    gsl_complex mesonPropagator = nonDiagonalMesonPropagator(
        solution.getParametersNJL(), 
        mottTemp, 
        solution.getUpQuarkChemicalPotential(), 
        solution.getDownQuarkChemicalPotential(), 
        solution.getStrangeQuarkChemicalPotential(), 
        mU, 
        mD, 
        mS, 
        k0, 
        0.0, 
        0.0, 
        mesonPropagatorPrecision, 
        mesonID
    );

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

void SU3NJL3DCutoffFixedChemPotTemp::findNondiagonalMesonMottTemperature(
    mesonState mesonIDAux,
    double precision, 
    MultiRootFindingMethod method, 
    double upQuarkEffectiveMassGuess, 
    double downQuarkEffectiveMassGuess, 
    double strangeQuarkEffectiveMassGuess,
    double mottTemperatureGuess
)
{   
    double x[4];
    x[0] = upQuarkEffectiveMassGuess; 
    x[1] = downQuarkEffectiveMassGuess; 
    x[2] = strangeQuarkEffectiveMassGuess;
    x[3] = mottTemperatureGuess;
    
    setMesonID(mesonIDAux);
    multiDimensionalRootFind(
        4, 
        precision, 
        &x[0], 
        this, 
        &SU3NJL3DCutoffNondiagonalMesonMottTemperatureFixedChemicalPotentials, 
        method
    );

    setUpQuarkEffectiveMass(x[0]);
    setDownQuarkEffectiveMass(x[1]);
    setStrangeQuarkEffectiveMass(x[2]);
    setTemperature(x[3]);
}

SU3NJL3DCutoffFixedChemPotTemp nondiagonalMesonMeltingPoint(
    SU3NJL3DCutoffVacuum vacuumSolution, 
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSolution, 
    mesonState mesonID, 
    double mesonPropertiesPrecision, 
    MultiRootFindingMethod method, 
    double mesonMassVacuumGuess, 
    double mesonWidthVacuumGuess
)
{   
    //calculate the meson mass behaviour for the provided range of NJL solutions
    vector<SU3NJL3DCutoffMeson> mesonFiniteT = mesonPropertiesFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
        vacuumSolution, 
        finiteTempSolution, 
        mesonID, 
        mesonPropertiesPrecision, 
        method, 
        mesonMassVacuumGuess, 
        mesonWidthVacuumGuess
    );

    //calculate Mott temperature of a given meson if it exists
    int meltingPointGuess = 0;
    double meltTest = 0;
    for (int i = 0; i < int(mesonFiniteT.size()); ++i)
    {   
        double k0AtMeltingPoint = mesonStateMassAtMeltingPoint(
            mesonFiniteT[i].getUpQuarkEffectiveMass(), 
            mesonFiniteT[i].getDownQuarkEffectiveMass(), 
            mesonFiniteT[i].getStrangeQuarkEffectiveMass(), 
            mesonID
        );
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
        mottSolution.findNondiagonalMesonMottTemperature(
            mesonID, 
            mesonPropertiesPrecision, 
            HYBRIDS, 
            mesonFiniteT[meltingPointGuess].getUpQuarkEffectiveMass(), 
            mesonFiniteT[meltingPointGuess].getDownQuarkEffectiveMass(), 
            mesonFiniteT[meltingPointGuess].getStrangeQuarkEffectiveMass(), 
            mesonFiniteT[meltingPointGuess].getTemperature()
        );
        cout << mottSolution.getTemperature() << "\n";
    }

    return mottSolution;
}

void SU3NJL3DCutoffFixedChemPotTemp::evaluateIsospinSymmetricCrossSections(
    SU3NJL3DCutoffParameters& parameters, 
    double precisionVacuum, 
    MultiRootFindingMethod methodVacuum, 
    double lightQuarkMassGuess, 
    double strangeQuarkMassGuess, 
    double nearVacuumTemperature, 
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
        nearVacuumTemperature, 
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

double SU3NJL3DCutoffFixedChemPotTemp::calculatePressure(double vacuumPressure)
{
    //Consider chemical potentials equal to effective chemical potentials
    //This holds if no vector interactions are considered
    //If vector interactions are considered, abort until this is updated
    LagrangianInteractions lagrangianInteractions = getParametersNJL().getDimensionfulCouplings().getLagrangianInteractions();
    if ( lagrangianInteractions!=SP4Q_DET2NFQ && lagrangianInteractions!=SP4Q_DET2NFQ_SP8Q  )
    {   
        cout << "Lagrangian interactions contain vector degrees of freedom! "
             << "The class SU3NJL3DCutoffFixedChemPotTemp is not prepared for this!\n";
        abort();
    }
    double upQuarkEffectiveChemicalPotential = upQuarkChemicalPotential;
    double downQuarkEffectiveChemicalPotential = downQuarkChemicalPotential;
    double strangeQuarkEffectiveChemicalPotential = strangeQuarkChemicalPotential;

    double pressureNJL = SU3NJL3DCutoffPressure(
        parametersNJL, 
        temperature,                                        
        upQuarkEffectiveMass,                                         
        downQuarkEffectiveMass,                                         
        strangeQuarkEffectiveMass,                                         
        upQuarkEffectiveChemicalPotential,                                         
        downQuarkEffectiveChemicalPotential,                                         
        strangeQuarkEffectiveChemicalPotential
    );

    pressureNJL = pressureNJL - vacuumPressure;

    return pressureNJL;
}

double SU3NJL3DCutoffFixedChemPotTemp::calculateEnergyDensity(double vacuumEnergyDensity)
{
    //Consider chemical potentials equal to effective chemical potentials
    //This holds if no vector interactions are considered
    //If vector interactions are considered, abort until this is updated
    LagrangianInteractions lagrangianInteractions = getParametersNJL().getDimensionfulCouplings().getLagrangianInteractions();
    if ( lagrangianInteractions!=SP4Q_DET2NFQ && lagrangianInteractions!=SP4Q_DET2NFQ_SP8Q  )
    {   
        cout << "Lagrangian interactions contain vector degrees of freedom! "
             << "The class SU3NJL3DCutoffFixedChemPotTemp is not prepared for this!\n";
        abort();
    }
    double upQuarkEffectiveChemicalPotential = upQuarkChemicalPotential;
    double downQuarkEffectiveChemicalPotential = downQuarkChemicalPotential;
    double strangeQuarkEffectiveChemicalPotential = strangeQuarkChemicalPotential;

    double energyNJL = SU3NJL3DCutoffEnergyDensity(
        parametersNJL, 
        temperature,                                           
        upQuarkEffectiveMass,                                            
        downQuarkEffectiveMass,                                            
        strangeQuarkEffectiveMass,                                            
        upQuarkEffectiveChemicalPotential,                                            
        downQuarkEffectiveChemicalPotential,                                            
        strangeQuarkEffectiveChemicalPotential
    );

    energyNJL = energyNJL - vacuumEnergyDensity;

    return energyNJL;
}


double SU3NJL3DCutoffFixedChemPotTemp::calculateEntropyDensity()
{
    //Consider chemical potentials equal to effective chemical potentials
    //This holds if no vector interactions are considered
    //If vector interactions are considered, abort until this is updated
    LagrangianInteractions lagrangianInteractions = getParametersNJL().getDimensionfulCouplings().getLagrangianInteractions();
    if ( lagrangianInteractions!=SP4Q_DET2NFQ && lagrangianInteractions!=SP4Q_DET2NFQ_SP8Q  )
    {   
        cout << "Lagrangian interactions contain vector degrees of freedom! "
             << "The class SU3NJL3DCutoffFixedChemPotTemp is not prepared for this!\n";
        abort();
    }
    double upQuarkEffectiveChemicalPotential = upQuarkChemicalPotential;
    double downQuarkEffectiveChemicalPotential = downQuarkChemicalPotential;
    double strangeQuarkEffectiveChemicalPotential = strangeQuarkChemicalPotential;

    double entropyNJL = SU3NJL3DCutoffEntropyDensity(
        parametersNJL, 
        temperature,                                             
        upQuarkEffectiveMass,                                              
        downQuarkEffectiveMass,                                              
        strangeQuarkEffectiveMass,                                              
        upQuarkEffectiveChemicalPotential,                                              
        downQuarkEffectiveChemicalPotential,                                              
        strangeQuarkEffectiveChemicalPotential
    );

    return entropyNJL;
}

void writeSolutionsToFile(vector<SU3NJL3DCutoffFixedChemPotTemp> solutions, string filename)
{
    int dataPrecision = 15;
    int colW = 25;

    std::ofstream file;
    file.open(filename, std::fstream::in | std::ofstream::out | std::ios::trunc);
    file.precision(dataPrecision);

    file.width(colW); file << "temperature[GeV]";

    file.width(colW); file << "effMassU[GeV]";
    file.width(colW); file << "effMassD[GeV]";
    file.width(colW); file << "effMassS[GeV]";

    file.width(colW); file << "chemPotU[GeV]";
    file.width(colW); file << "chemPotD[GeV]";
    file.width(colW); file << "chemPotS[GeV]";

    file.width(colW); file << "pressure[GeV^4]";
    file.width(colW); file << "energyDensity[GeV^4]";
    file.width(colW); file << "entropyDensity[GeV^3]";
    file << "\n";

    for (int i = 0; i < int(solutions.size()); ++i)
    {   
        file.width(25);   file << solutions[i].getTemperature();
        
        file.width(25);   file << solutions[i].getUpQuarkEffectiveMass();
        file.width(25);   file << solutions[i].getDownQuarkEffectiveMass();
        file.width(25);   file << solutions[i].getStrangeQuarkEffectiveMass();

        file.width(25);   file << solutions[i].getUpQuarkChemicalPotential();
        file.width(25);   file << solutions[i].getDownQuarkChemicalPotential();
        file.width(25);   file << solutions[i].getStrangeQuarkChemicalPotential();

        file.width(25);   file << solutions[i].getPressure();
        file.width(25);   file << solutions[i].getEnergyDensity();
        file.width(25);   file << solutions[i].getEntropyDensity();

        file << "\n";
    }

    file.close();
}

void SU3NJL3DCutoffFixedChemPotTemp::evaluateInMediumMassesAndThermodynamics(
	SU3NJL3DCutoffParameters& parameters, 
    double precisionVacuum, 
    MultiRootFindingMethod methodVacuum, 
    double upQuarkMassGuess, 
    double downQuarkMassGuess, 
    double strangeQuarkMassGuess, 
    double nearVacuumTemperature,
	double temperature, 
	int numberOfPoints, 
    double precisionVacToFinTemp, 
    MultiRootFindingMethod methodVacToFinTemp
)
{
	SU3NJL3DCutoffVacuum vacuum = SU3NJL3DCutoffVacuum::calculateVacuumMasses(
        parameters,                                    
        precisionVacuum,                                    
        methodVacuum,                                    
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess
    );

    //solve model at zero chemical potential up to some finite temperature
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTSolution = solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
        vacuum, 
        nearVacuumTemperature,
        temperature, 
        numberOfPoints, 
        precisionVacToFinTemp, 
        methodVacToFinTemp
    );

    //calculate vacuum pressure
    cout << "\nCalculating the vacuum pressure...\n";
    double pressureVac = vacuum.calculatePressure();
    //cout << "vacuumPressure[GeV^4]=" << pressureVac << "\n";

    cout << "\nCalculating thermodynamics for the found solutions...\n";
    for (int i = 0; i < int(finiteTSolution.size()); ++i)
    {   
        finiteTSolution[i].setPressure(pressureVac);
        finiteTSolution[i].setEnergyDensity(-pressureVac);
        finiteTSolution[i].setEntropyDensity();
    }

    //Add solution with vacuum solution to start of the vector if near vacuum temperature is not zero
    if (nearVacuumTemperature>0)
    {
        SU3NJL3DCutoffFixedChemPotTemp auxVacuum(vacuum);
        finiteTSolution.insert(finiteTSolution.begin(), auxVacuum);
    }

    string filename = "SU3NJL3DCutoffFixedChemPotTemp";
    filename = filename + "_" + finiteTSolution[0].getParametersNJL().getParameterSetName();
    filename = filename + "_TMin" + to_string(finiteTSolution[0].getTemperature());
    filename = filename + "_TMax" + to_string(finiteTSolution[finiteTSolution.size()-1].getTemperature());
    filename = filename + "_CP0";
    replaceChar(filename, '.', 'p');
    filename =  filename +".dat";
    writeSolutionsToFile(finiteTSolution, filename);
}

void calculateThermodynamics(
    SU3NJL3DCutoffVacuum& vacuum, 
    vector<SU3NJL3DCutoffFixedChemPotTemp>& solution
)
{
    //calculate vacuum pressure
    cout << "\nCalculating the vacuum pressure...\n";
    double vacuumPressure = vacuum.calculatePressure();
    //cout << "vacuumPressure[GeV^4]=" << vacuumPressure << "\n";

    cout << "\nCalculating thermodynamics for the provided solutions...\n";
    for (int i = 0; i < int(solution.size()); ++i)
    {   
        solution[i].setPressure(vacuumPressure);
        solution[i].setEnergyDensity(-vacuumPressure);
        solution[i].setEntropyDensity();
    }
}

void SU3NJL3DCutoffFixedChemPotTemp::computeThermoFixedChemPotTrajectory(
    SU3NJL3DCutoffParameters& parameters, 
    double precisionVacuum,
    MultiRootFindingMethod methodVacuum,
    double upQuarkMassGuess,
    double downQuarkMassGuess,
    double strangeQuarkMassGuess,
    double chemicalPotential,
    int numberOfPointsVacToChemPot,
    double precisionVacToChemPot,
    MultiRootFindingMethod methodVacToChemPot,
    double temperature,
    int numberOfPointsUpToTemp,
	double precisionUpToTemp,
	MultiRootFindingMethod methodUpToTemp,
    string customSuffix
)
{
    SU3NJL3DCutoffVacuum vacuum = SU3NJL3DCutoffVacuum::calculateVacuumMasses(
        parameters,                                    
        precisionVacuum,                                    
        methodVacuum,                                    
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess
    );

    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteChemPotSolution = solveVacuumToFiniteChemicalPotential(
        vacuum, 
        chemicalPotential, 
        numberOfPointsVacToChemPot, 
        precisionVacToChemPot, 
        methodVacToChemPot
    );

    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTempSolution = solveUpToTemperature(
        finiteChemPotSolution[finiteChemPotSolution.size()-1], 
        temperature, 
        numberOfPointsUpToTemp, 
        precisionUpToTemp, 
        methodUpToTemp
    );

    calculateThermodynamics(vacuum, finiteTempSolution);

    // coutput to file
    string filename = "SU3NJL3DCutoffFixedChemPotTemp";
    if( !customSuffix.empty() )
    {
        filename = filename + "_" + customSuffix;
    }
    else
    {
        filename = filename + "_" + vacuum.getParametersNJL().getParameterSetName();
        filename = filename + "_TMin" + trim0ToDot0(finiteTempSolution[0].getTemperature());
        filename = filename + "_TMax" + trim0ToDot0(finiteTempSolution[finiteTempSolution.size()-1].getTemperature());
        filename = filename + "_CPU"  + trim0ToDot0(finiteTempSolution[0].getUpQuarkChemicalPotential());
    }
    replaceChar(filename, '.', 'p');
    filename =  filename + ".dat";
    writeSolutionsToFile(finiteTempSolution, filename);
}

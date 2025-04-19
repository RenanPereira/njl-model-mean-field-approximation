#include <cmath>
#include <iostream>
#include "physics_utils/physical_constants.h"
#include "njl_model/n_fermion_line_integrals/one_fermion_line_integral_3d_cutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedTempRhoBEqualChemPot.h"
#include "utils/format_utils.h"

using namespace std;


SU3NJL3DCutoffFixedTempRhoBEqualChemPot::SU3NJL3DCutoffFixedTempRhoBEqualChemPot(void* auxiliar)
{	
	parametersNJL = ((class SU3NJL3DCutoffFixedTempRhoBEqualChemPot *)(auxiliar))->parametersNJL;

	upQuarkEffectiveMass = ((class SU3NJL3DCutoffFixedTempRhoBEqualChemPot *)(auxiliar))->upQuarkEffectiveMass;
	downQuarkEffectiveMass = ((class SU3NJL3DCutoffFixedTempRhoBEqualChemPot *)(auxiliar))->downQuarkEffectiveMass;
	strangeQuarkEffectiveMass = ((class SU3NJL3DCutoffFixedTempRhoBEqualChemPot *)(auxiliar))->strangeQuarkEffectiveMass;
	quarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffFixedTempRhoBEqualChemPot *)(auxiliar))->quarkEffectiveChemicalPotential;

	temperature = ((class SU3NJL3DCutoffFixedTempRhoBEqualChemPot *)(auxiliar))->temperature;
	baryonDensity = ((class SU3NJL3DCutoffFixedTempRhoBEqualChemPot *)(auxiliar))->baryonDensity;
};


SU3NJL3DCutoffFixedTempRhoBEqualChemPot::SU3NJL3DCutoffFixedTempRhoBEqualChemPot(SU3NJL3DCutoffParameters parametersNJLAux)
{
	parametersNJL = parametersNJLAux;
}


SU3NJL3DCutoffFixedTempRhoBEqualChemPot::SU3NJL3DCutoffFixedTempRhoBEqualChemPot(
    SU3NJL3DCutoffParameters parametersNJLAux,                                                                                 
    double temperatureAux
)
{
	parametersNJL = parametersNJLAux;
    temperature = temperatureAux;
}


SU3NJL3DCutoffFixedTempRhoBEqualChemPot::SU3NJL3DCutoffFixedTempRhoBEqualChemPot(
    SU3NJL3DCutoffParameters parametersNJLAux,                                                                             
    double temperatureAux,                                                                             
    double baryonDensityAux
)
{
	parametersNJL = parametersNJLAux;
    temperature = temperatureAux;
    baryonDensity = baryonDensityAux;
}


SU3NJL3DCutoffFixedTempRhoBEqualChemPot::SU3NJL3DCutoffFixedTempRhoBEqualChemPot(
    SU3NJL3DCutoffParameters parametersNJLAux,         
    double temperatureAux,                                                                             
    double baryonDensityAux,                                                                             
    double upQuarkEffectiveMassAux,                                                                             
    double downQuarkEffectiveMassAux,                                                                             
    double strangeQuarkEffectiveMassAux,                                                                             
    double quarkEffectiveChemicalPotentialAux
)
{
	parametersNJL = parametersNJLAux;

    temperature = temperatureAux;
    baryonDensity = baryonDensityAux;
    
    upQuarkEffectiveMass = upQuarkEffectiveMassAux;
    downQuarkEffectiveMass = downQuarkEffectiveMassAux;
    strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux;
    quarkEffectiveChemicalPotential = quarkEffectiveChemicalPotentialAux;
}


SU3NJL3DCutoffFixedTempRhoBEqualChemPot::SU3NJL3DCutoffFixedTempRhoBEqualChemPot(SU3NJL3DCutoffVacuum vacuum)
{   
    parametersNJL = vacuum.getParametersNJL();

    temperature = 0.0;
    baryonDensity = 0.0;

    upQuarkEffectiveMass = vacuum.getUpQuarkEffectiveMass();
    downQuarkEffectiveMass = vacuum.getDownQuarkEffectiveMass();
    strangeQuarkEffectiveMass = vacuum.getStrangeQuarkEffectiveMass();
    quarkEffectiveChemicalPotential = 0.0;

    upQuarkChemicalPotential = 0.0;
	downQuarkChemicalPotential = 0.0;
	strangeQuarkChemicalPotential = 0.0;
	baryonChemicalPotential = 0.0;

    setSigmasDensitiesChemicalPotentials(
        upQuarkEffectiveMass,                                  
        downQuarkEffectiveMass,                                  
        strangeQuarkEffectiveMass,                                  
        quarkEffectiveChemicalPotential,                                  
        quarkEffectiveChemicalPotential,                                  
        quarkEffectiveChemicalPotential
    );

    pressure = 0.0;
    energyDensity = 0.0;
    entropyDensity = 0.0;
}


void SU3NJL3DCutoffFixedTempRhoBEqualChemPot::solve(double precision, MultiRootFindingMethod method, 
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
    
    multiDimensionalRootFind(4, precision, &x[0], this, &gapEquations, method);

	setUpQuarkEffectiveMass(x[0]);
	setDownQuarkEffectiveMass(x[1]);
	setStrangeQuarkEffectiveMass(x[2]);
	setQuarkEffectiveChemicalPotential(x[3]);
}


int SU3NJL3DCutoffFixedTempRhoBEqualChemPot::gapEquations(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{
    if (x->size != 4) 
    {   
        string functionName = "SU3NJL3DCutoffFixedTempRhoBEqualChemPot::gapEquations";
        cout << "Error: gsl_vector 'x' has insufficient size in " << functionName << "\n";
        cout << "(expected 4, got " << x->size << ").\n";
        return GSL_FAILURE;
    }

	//Define variables
    double mU = gsl_vector_get(x,0);
    double mD = gsl_vector_get(x,1);
    double mS = gsl_vector_get(x,2);
    double effCP = gsl_vector_get(x,3);

    //Get parameters
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot solution(auxiliar);
    NJLDimensionfulCouplings couplings = solution.getParametersNJL().getDimensionfulCouplings();
    double m0U = solution.getParametersNJL().getUpQuarkCurrentMass();
    double m0D = solution.getParametersNJL().getDownQuarkCurrentMass();
    double m0S = solution.getParametersNJL().getStrangeQuarkCurrentMass();
    double rhoB = solution.getBaryonDensity();

    //Calculate the sigma field and density for each quark flavour
    solution.setSigmasDensitiesChemicalPotentials(mU, mD, mS, effCP, effCP, effCP);
    double sigmaU = solution.getUpQuarkSigma();
    double sigmaD = solution.getDownQuarkSigma();
    double sigmaS = solution.getStrangeQuarkSigma();
    double rhoU = solution.getUpQuarkDensity();
    double rhoD = solution.getDownQuarkDensity();
    double rhoS = solution.getStrangeQuarkDensity();

	//System of equations
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


bool SU3NJL3DCutoffFixedTempRhoBEqualChemPot::testSolution(double precision)
{   
    double x[4];
    x[0] = getUpQuarkEffectiveMass();
    x[1] = getDownQuarkEffectiveMass();
    x[2] = getStrangeQuarkEffectiveMass();
    x[3] = getQuarkEffectiveChemicalPotential();

    //The test below (gsl) resturns 0 if the sum_i abs(residual_i) < precision
    int gslTest = multiDimensionalRootFindTestResidual(4, precision, &x[0], this, &gapEquations);

    if (gslTest==0){ return true; }
    else{ return false; }
}


double SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculatePressure(double vacuumPressure)
{
    double pressureNJL = SU3NJL3DCutoffPressure(
        parametersNJL, temperature,                                        
        upQuarkEffectiveMass,                                         
        downQuarkEffectiveMass,                                         
        strangeQuarkEffectiveMass,                                         
        quarkEffectiveChemicalPotential,                                         
        quarkEffectiveChemicalPotential,                                         
        quarkEffectiveChemicalPotential
    );

    pressureNJL = pressureNJL - vacuumPressure;

    return pressureNJL;
}


double SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateEnergyDensity(double vacuumEnergyDensity)
{
    double energyNJL = SU3NJL3DCutoffEnergyDensity(
        parametersNJL, 
        temperature,                                           
        upQuarkEffectiveMass,                                            
        downQuarkEffectiveMass,                                            
        strangeQuarkEffectiveMass,                                            
        quarkEffectiveChemicalPotential,                                            
        quarkEffectiveChemicalPotential,                                            
        quarkEffectiveChemicalPotential
    );

    energyNJL = energyNJL - vacuumEnergyDensity;

    return energyNJL;
}


double SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateEntropyDensity()
{
    double entropyNJL = SU3NJL3DCutoffEntropyDensity(
        parametersNJL, temperature,                                             
        upQuarkEffectiveMass,                                              
        downQuarkEffectiveMass,                                              
        strangeQuarkEffectiveMass,                                              
        quarkEffectiveChemicalPotential,                                              
        quarkEffectiveChemicalPotential,                                              
        quarkEffectiveChemicalPotential
    );

    return entropyNJL;
}


void SU3NJL3DCutoffFixedTempRhoBEqualChemPot::setSigmasDensitiesChemicalPotentials(
    double effMassU,
    double effMassD,
    double effMassS,
    double effCPU,
    double effCPD,
    double effCPS
)
{
    //Calculate the outputs of this solution
    NJLDimensionfulCouplings couplings = parametersNJL.getDimensionfulCouplings();

    NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();

    double Nc = parametersNJL.getNumberOfColours();
    
    double sigmaIntegralPrecision = parametersNJL.getSigmaIntegralPrecision();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

    double T = getTemperature();
   
    //Calculate the sigma field and density for each quark flavour
    double sigmaU = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCPU, effMassU, sigmaIntegralPrecision);
    double sigmaD = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCPD, effMassD, sigmaIntegralPrecision);
    double sigmaS = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCPS, effMassS, sigmaIntegralPrecision);

    double rhoU = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCPU, effMassU, thermoIntegralPrecision);
    double rhoD = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCPD, effMassD, thermoIntegralPrecision);
    double rhoS = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCPS, effMassS, thermoIntegralPrecision);

    //Calculate quark chemical potential from the effective ones
    double cPU = SU3NJLQuarkChemicalPotential(couplings, effCPU, rhoU, rhoD, rhoS, sigmaU, sigmaD, sigmaS);
    double cPD = SU3NJLQuarkChemicalPotential(couplings, effCPD, rhoD, rhoS, rhoU, sigmaD, sigmaS, sigmaU);
    double cPS = SU3NJLQuarkChemicalPotential(couplings, effCPS, rhoS, rhoU, rhoD, sigmaS, sigmaU, sigmaD);

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
}


void SU3NJL3DCutoffFixedTempRhoBEqualChemPot::setSigmasDensitiesChemicalPotentials()
{    
    double effMassU = getUpQuarkEffectiveMass();
    double effMassD = getDownQuarkEffectiveMass();
    double effMassS = getStrangeQuarkEffectiveMass();
    double effCP = getQuarkEffectiveChemicalPotential();

    setSigmasDensitiesChemicalPotentials(
        effMassU,
        effMassD,
        effMassS,
        effCP,
        effCP,
        effCP
    );
}


vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> SU3NJL3DCutoffFixedTempRhoBEqualChemPot::solveFromVacuumToFiniteBaryonDensity(
    SU3NJL3DCutoffVacuum vacuum, 
    double minimumBaryonDensity, double maximumBaryonDensity, int numberOfPoints, 
    double gapPrecision, MultiRootFindingMethod method, bool storeToFile
)
{   
    //Analyse vacuum solution
    double pressureVacuum = vacuum.calculatePressure();
    double energyVacuum = vacuum.calculateEnergyDensity();
    cout << "pressureVacuum=" << pressureVacuum << "\n";
    cout << "energyDensityVacuum=" << energyVacuum << "\n";
    
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot vacuumSol(vacuum);

    //Build guesses from the vacuum solution
    double mUGuess = vacuumSol.getUpQuarkEffectiveMass();
    double mDGuess = vacuumSol.getDownQuarkEffectiveMass();
    double mSGuess = vacuumSol.getStrangeQuarkEffectiveMass();
    double effCPGuess = vacuumSol.getUpQuarkEffectiveMass()*(1+1E-4);

    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> solutions;
    solutions.push_back(vacuumSol);
    double delta = (maximumBaryonDensity-minimumBaryonDensity)/(numberOfPoints-1);
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        //Step in density
        double rhoB = minimumBaryonDensity + i*delta;

        SU3NJL3DCutoffFixedTempRhoBEqualChemPot inMediumSol(
            vacuum.getParametersNJL(),
            0.0,
            rhoB
        );

        //Find quark masses and effective chemical potential
        inMediumSol.solve(gapPrecision, method, mUGuess, mDGuess, mSGuess, effCPGuess);

        //Guesses for next step
        mUGuess = inMediumSol.getUpQuarkEffectiveMass();
        mDGuess = inMediumSol.getDownQuarkEffectiveMass();
        mSGuess = inMediumSol.getStrangeQuarkEffectiveMass();
        effCPGuess = inMediumSol.getQuarkEffectiveChemicalPotential();  

        //Calculate the sigma field and density for each quark flavour
        inMediumSol.setSigmasDensitiesChemicalPotentials();
     
        //Calculate thermodynamics
        inMediumSol.setPressure( pressureVacuum );
        inMediumSol.setEnergyDensity( energyVacuum );
        inMediumSol.setEntropyDensity();

        //Push to solutions vector if solution passes the test
        bool testSolution = inMediumSol.testSolution(gapPrecision);
        if (testSolution==true){ solutions.push_back(inMediumSol); }
        else{ cout << "Solution fails test!"; }

        //Print solution to console
        printf("rhoB=%.6f, mU=%.6f, mD=%.6f, mS=%.6f, effCP=%.6f [GeV] \n", 
                rhoB/pow(hc_GeVfm, 3),
                inMediumSol.getUpQuarkEffectiveMass(), 
                inMediumSol.getDownQuarkEffectiveMass(), 
                inMediumSol.getStrangeQuarkEffectiveMass(),
                inMediumSol.getQuarkEffectiveChemicalPotential());
    }
    
    if ( storeToFile )
    {
        //Store calculation to file
        string filename = string("SU3NJL3DCutoffEqualChemPot_") + vacuum.getParametersNJL().getParameterSetName()
                                                                + "T0.0"
                                                                + "rhoBMin" + trim0ToDot0(minimumBaryonDensity/pow(hc_GeVfm,3))
                                                                + "rhoBMax" + trim0ToDot0(maximumBaryonDensity/pow(hc_GeVfm,3))
                                                                + "N" + to_string(numberOfPoints)
                                                                + ".dat";
        
        SU3NJL3DCutoffFixedTempRhoBEqualChemPot::writeToFile(solutions, filename, true);
    }

    return solutions;
}


void SU3NJL3DCutoffFixedTempRhoBEqualChemPot::writeToFile(vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> solutions, string fileName, bool columnsDescription)
{
    //Create file to store solution
    std::ofstream aux_file_NJL;
    aux_file_NJL.open(fileName, std::ofstream::out | std::ios::trunc);
    aux_file_NJL.precision(15);

    if ( columnsDescription )
    {
        aux_file_NJL.width(25);   aux_file_NJL << "rhoB[fm-3]";
        aux_file_NJL.width(25);   aux_file_NJL << "T[GeV]";
        aux_file_NJL.width(25);   aux_file_NJL << "massU[GeV]";
        aux_file_NJL.width(25);   aux_file_NJL << "massD[GeV]";
        aux_file_NJL.width(25);   aux_file_NJL << "massS[GeV]";
        aux_file_NJL.width(25);   aux_file_NJL << "quarkChemPot[GeV]";
        aux_file_NJL.width(25);   aux_file_NJL << "baryonChemPot[GeV]";
        aux_file_NJL.width(25);   aux_file_NJL << "pressure[GeVfm-3]";
        aux_file_NJL.width(25);   aux_file_NJL << "energyDensity[GeVfm-3]";
        aux_file_NJL.width(25);   aux_file_NJL << "entropyDensity[fm-3]";
        aux_file_NJL.width(25);   aux_file_NJL << "sigmaU[GeV]";
        aux_file_NJL.width(25);   aux_file_NJL << "sigmaD[GeV]";
        aux_file_NJL.width(25);   aux_file_NJL << "sigmaS[GeV]";
        aux_file_NJL.width(25);   aux_file_NJL << "rhoU[fm-3]";
        aux_file_NJL.width(25);   aux_file_NJL << "rhoD[fm-3]";
        aux_file_NJL.width(25);   aux_file_NJL << "rhoS[fm-3]";
        aux_file_NJL << std::endl;
    }

    for (int i = 0; i < int( solutions.size() ); ++i)
    {   
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getBaryonDensity()/pow(hc_GeVfm,3);
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getTemperature();
   		aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getUpQuarkEffectiveMass();
   		aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getDownQuarkEffectiveMass();
   		aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getStrangeQuarkEffectiveMass();
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getQuarkEffectiveChemicalPotential();
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getBaryonChemicalPotential();
   		aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getPressure()/pow(hc_GeVfm,3);
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getEnergyDensity()/pow(hc_GeVfm,3);
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getEntropyDensity()/pow(hc_GeVfm,3);
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getUpQuarkSigma();
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getDownQuarkSigma();
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getStrangeQuarkSigma();
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getUpQuarkDensity();
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getDownQuarkDensity();
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getStrangeQuarkDensity();
		aux_file_NJL << std::endl;
    }
}


SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint::ChiralTransitionPoint(
	SU3NJL3DCutoffParameters params,
	double temp,
	double mU_broken, double mD_broken, double mS_broken, double effCP_broken,
	double mU_restored, double mD_restored, double mS_restored, double effCP_restored
)
	: upQuarkEffectiveMassBroken(mU_broken),
	  downQuarkEffectiveMassBroken(mD_broken),
	  strangeQuarkEffectiveMassBroken(mS_broken),
	  quarkEffectiveChemicalPotentialBroken(effCP_broken),
	  upQuarkEffectiveMassRestored(mU_restored),
	  downQuarkEffectiveMassRestored(mD_restored),
	  strangeQuarkEffectiveMassRestored(mS_restored),
	  quarkEffectiveChemicalPotentialRestored(effCP_restored),
	  parametersNJL(params),
	  temperature(temp)
{}


int SU3NJL3DCutoffFixedTempRhoBEqualChemPot::chiralTransitionEquations(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{
    if (x->size != 8) 
    {   
        string functionName = "SU3NJL3DCutoffFixedTempRhoBEqualChemPot::chiralTransitionEquations";
        cout << "Error: gsl_vector 'x' has insufficient size in " << functionName << "\n";
        cout << "(expected 4, got " << x->size << ").\n";
        return GSL_FAILURE;
    }

    //Define variables
    double mU_broken = gsl_vector_get(x,0);
    double mD_broken = gsl_vector_get(x,1);
    double mS_broken = gsl_vector_get(x,2);
    double effCP_broken = gsl_vector_get(x,3);

    double mU_restored = gsl_vector_get(x,4);
    double mD_restored = gsl_vector_get(x,5);
    double mS_restored = gsl_vector_get(x,6);
    double effCP_restored = gsl_vector_get(x,7);

    //Read parameters from source object
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot source(auxiliar);

    NJLDimensionfulCouplings couplings = source.getParametersNJL().getDimensionfulCouplings();
    double m0U = source.getParametersNJL().getUpQuarkCurrentMass();
    double m0D = source.getParametersNJL().getDownQuarkCurrentMass();
    double m0S = source.getParametersNJL().getStrangeQuarkCurrentMass();
    double T = source.getTemperature();

    //Broken phase
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot brokenPhase(
        source.getParametersNJL(),
        T,
        0.0,
        mU_broken,
        mD_broken,
        mS_broken,
        effCP_broken
    );
    brokenPhase.setSigmasDensitiesChemicalPotentials();

    double sigmaU_broken = brokenPhase.getUpQuarkSigma();
    double sigmaD_broken = brokenPhase.getDownQuarkSigma();
    double sigmaS_broken = brokenPhase.getStrangeQuarkSigma();

    //Get all the necessary quantities
    double rhoU_broken = brokenPhase.getUpQuarkDensity();
    double rhoD_broken = brokenPhase.getDownQuarkDensity();
    double rhoS_broken = brokenPhase.getStrangeQuarkDensity();

    double baryonCP_broken = brokenPhase.getBaryonChemicalPotential();
    double pressure_broken = brokenPhase.calculatePressure(0.0);//the vacuum values can be any value

    //Restored phase
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot restoredPhase(
        source.getParametersNJL(),
        T,
        0.0,
        mU_restored,
        mD_restored,
        mS_restored,
        effCP_restored
    );
    restoredPhase.setSigmasDensitiesChemicalPotentials();

    double sigmaU_restored = restoredPhase.getUpQuarkSigma();
    double sigmaD_restored = restoredPhase.getDownQuarkSigma();
    double sigmaS_restored = restoredPhase.getStrangeQuarkSigma();

    //Get all the necessary quantities
    double rhoU_restored = restoredPhase.getUpQuarkDensity();
    double rhoD_restored = restoredPhase.getDownQuarkDensity();
    double rhoS_restored = restoredPhase.getStrangeQuarkDensity();

    double baryonCP_restored = restoredPhase.getBaryonChemicalPotential();
    double pressure_restored = restoredPhase.calculatePressure(0.0);//the vacuum values can be any value

    //Equations in the broken phase
    double f0 = SU3NJLNulledGapEquation(couplings, mU_broken-m0U, sigmaU_broken, sigmaD_broken, sigmaS_broken, rhoU_broken, rhoD_broken, rhoS_broken);
    double f1 = SU3NJLNulledGapEquation(couplings, mD_broken-m0D, sigmaD_broken, sigmaS_broken, sigmaU_broken, rhoD_broken, rhoS_broken, rhoU_broken);
    double f2 = SU3NJLNulledGapEquation(couplings, mS_broken-m0S, sigmaS_broken, sigmaU_broken, sigmaD_broken, rhoS_broken, rhoU_broken, rhoD_broken);

    //Equations in the restored phase
    double f3 = SU3NJLNulledGapEquation(couplings, mU_restored-m0U, sigmaU_restored, sigmaD_restored, sigmaS_restored, rhoU_restored, rhoD_restored, rhoS_restored);
    double f4 = SU3NJLNulledGapEquation(couplings, mD_restored-m0D, sigmaD_restored, sigmaS_restored, sigmaU_restored, rhoD_restored, rhoS_restored, rhoU_restored);
    double f5 = SU3NJLNulledGapEquation(couplings, mS_restored-m0S, sigmaS_restored, sigmaU_restored, sigmaD_restored, rhoS_restored, rhoU_restored, rhoD_restored);
    
    //Maxwell construction
    double f6 = baryonCP_broken - baryonCP_restored;
    double f7 = pressure_broken - pressure_restored;

    gsl_vector_set (f,0,f0);
    gsl_vector_set (f,1,f1);
    gsl_vector_set (f,2,f2);
    gsl_vector_set (f,3,f3);
    gsl_vector_set (f,4,f4);
    gsl_vector_set (f,5,f5);
    gsl_vector_set (f,6,f6);
    gsl_vector_set (f,7,f7);

    return GSL_SUCCESS;
}


SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateChiralTransitionPoint(
    double T,
    const ChiralTransitionPoint& guess,
    double precision,
    MultiRootFindingMethod method) 
{
    double x[8] = {
        guess.upQuarkEffectiveMassBroken,
        guess.downQuarkEffectiveMassBroken,
        guess.strangeQuarkEffectiveMassBroken,
        guess.quarkEffectiveChemicalPotentialBroken,
        guess.upQuarkEffectiveMassRestored,
        guess.downQuarkEffectiveMassRestored,
        guess.strangeQuarkEffectiveMassRestored,
        guess.quarkEffectiveChemicalPotentialRestored
    };

    SU3NJL3DCutoffFixedTempRhoBEqualChemPot aux(guess.parametersNJL, T);
    multiDimensionalRootFind(8, precision, x, &aux, &chiralTransitionEquations, method);

    return ChiralTransitionPoint (
        guess.parametersNJL, T,
        x[0], x[1], x[2], x[3], 
        x[4], x[5], x[6], x[7] 
    );
}


vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint> SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateFirstOrderLine(
    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> solutions, 
    double precision, 
    MultiRootFindingMethod method,
    double deltaT,
    double massDifferenceCEP
)
{       
    //Run check that solutions have the same temperature
    for (int i = 0; i < int(solutions.size()); ++i)
    {   
        if ( fabs(solutions[1].getTemperature()-solutions[i].getTemperature())>0 ){ cout << "Temperature is not fixed in this vector of solutions!\n"; abort(); }
    }

    SU3NJL3DCutoffParameters parametersNJL = solutions[0].getParametersNJL();
    double T = solutions[0].getTemperature();
    
    //Finding guesses for transition point
    int first = 0;
    for (int i = 1; i < int(solutions.size()); ++i)
    {
        double delta = solutions[i-1].getQuarkEffectiveChemicalPotential() - solutions[i].getQuarkEffectiveChemicalPotential();
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
            double delta = solutions[i-1].getQuarkEffectiveChemicalPotential() - solutions[i].getQuarkEffectiveChemicalPotential();
            if ( delta<0 )
            {
                second = i;
                //cout << second << "\n";
                break;
            }
        }
    }

    vector<ChiralTransitionPoint> firtOrderLine;
    if ( first>0 && second>0 )
    {
        //Estimate guesses
        double mU_broken = solutions[first].getUpQuarkEffectiveMass();
        double mD_broken = solutions[first].getDownQuarkEffectiveMass();
        double mS_broken = solutions[first].getStrangeQuarkEffectiveMass();
        double effCP_broken = solutions[first].getQuarkEffectiveChemicalPotential();
        double mU_restored = solutions[second].getUpQuarkEffectiveMass();
        double mD_restored = solutions[second].getDownQuarkEffectiveMass();
        double mS_restored = solutions[second].getStrangeQuarkEffectiveMass();
        double effCP_restored = solutions[second].getQuarkEffectiveChemicalPotential();

        ChiralTransitionPoint point = {
            parametersNJL,
            T,
            mU_broken,
            mD_broken,
            mS_broken,
            0.5*effCP_broken + 0.5*effCP_restored,
            mU_restored,
            mD_restored,
            mS_restored,
            0.5*effCP_broken + 0.5*effCP_restored
        };
        
        point = calculateChiralTransitionPoint(T, point, precision, method);

        mU_broken = point.upQuarkEffectiveMassBroken; 
        mD_broken = point.downQuarkEffectiveMassBroken;
        mS_broken = point.strangeQuarkEffectiveMassBroken;
        effCP_broken = point.quarkEffectiveChemicalPotentialBroken;
        mU_restored = point.upQuarkEffectiveMassRestored; 
        mD_restored = point.downQuarkEffectiveMassRestored; 
        mS_restored = point.strangeQuarkEffectiveMassRestored;
        effCP_restored = point.quarkEffectiveChemicalPotentialRestored;
        
        //Store the chiral transition point
        firtOrderLine.push_back(point);
        
        printf("Temperature=%.6f [GeV]\n", T);
        printf("Broken phase:\n");
        printf("mU=%.6f, mD=%.6f, mS=%.6f, effCP=%.6f [GeV] \n", mU_broken, mD_broken, mS_broken, effCP_broken);
        printf("Restored phase:\n");
        printf("mU=%.6f, mD=%.6f, mS=%.6f, effCP=%.6f [GeV] \n\n", mU_restored, mD_restored, mS_restored, effCP_restored);

        //Calculate if there is a first order phase transition by checking if the difference between masses in both phase is large
		double deltaCEP = fabs(mU_broken - mU_restored);

        //Find first order line by increasing temperature while the mass difference is not reached: repeat the above calculation
        T = 0.001;//Set temperature to 1GeV for the first step
        while ( deltaCEP>massDifferenceCEP )
        {
            point = calculateChiralTransitionPoint(T, point, precision, method);

            mU_broken = point.upQuarkEffectiveMassBroken; 
            mD_broken = point.downQuarkEffectiveMassBroken;
            mS_broken = point.strangeQuarkEffectiveMassBroken;
            effCP_broken = point.quarkEffectiveChemicalPotentialBroken;
            mU_restored = point.upQuarkEffectiveMassRestored; 
            mD_restored = point.downQuarkEffectiveMassRestored; 
            mS_restored = point.strangeQuarkEffectiveMassRestored;
            effCP_restored = point.quarkEffectiveChemicalPotentialRestored;

            //Store the chiral transition point
            firtOrderLine.push_back(point);
            
            printf("Temperature=%.6f [GeV]\n", T);
            printf("Broken phase:\n");
            printf("mU=%.6f, mD=%.6f, mS=%.6f, effCP=%.6f [GeV] \n", mU_broken, mD_broken, mS_broken, effCP_broken);
            printf("Restored phase:\n");
            printf("mU=%.6f, mD=%.6f, mS=%.6f, effCP=%.6f [GeV] \n\n", mU_restored, mD_restored, mS_restored, effCP_restored);
            
            deltaCEP = abs(mU_broken - mU_restored);
            
            //Increase temperature
            T = T + deltaT;
        }
    }
    else
    {
        cout << "No first-order phase transition point was detected at zero temperature!\n";
    }

    return firtOrderLine;
}


vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint> SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateFirstOrderLine(
    SU3NJL3DCutoffVacuum vacuum,
    double minimumBaryonDensity,
    double maximumBaryonDensity,
    int numberOfPoints,
    double precisionZeroTempSol,
    MultiRootFindingMethod methodZeroTempSol,
    bool storeZeroTempSolToFile,
    double precisionTransitionPointSol, 
    MultiRootFindingMethod methodTransitionPointSol,
    double deltaT,
    double massDifferenceCEP
)
{    
    //Solve the model at zero temperature in order to find guess for chiral transition
    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> zeroTempSol = solveFromVacuumToFiniteBaryonDensity(
        vacuum, 
        minimumBaryonDensity, 
        maximumBaryonDensity, 
        numberOfPoints, 
        precisionZeroTempSol, 
        methodZeroTempSol, 
        storeZeroTempSolToFile
    );
    
    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint> firstOrderLine = 
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateFirstOrderLine(
        zeroTempSol, 
        precisionTransitionPointSol, 
        methodTransitionPointSol, 
        deltaT, 
        massDifferenceCEP
    );

    return firstOrderLine;
}


void SU3NJL3DCutoffFixedTempRhoBEqualChemPot::writeToFile(
    SU3NJL3DCutoffVacuum vacuum,
    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint> firstOrderLine,
    string fileName, 
    bool columnsDescription
)
{
    // Get Vacuum Pressure
    double pressureVacuum = vacuum.calculatePressure();

    //Create file to store first order line
    std::ofstream file;
    file.open(fileName, std::ofstream::out | std::ios::trunc);
    file.precision(15);

    if ( columnsDescription )
    {
        file.width(25);   file << "quarkChemPot[GeV]";
        file.width(25);   file << "baryonChemPot[GeV]";
        file.width(25);   file << "T[GeV]";
        file.width(25);   file << "pressure[GeVfm-3]";

        file.width(25);   file << "rhoB_broken[fm-3]";
        file.width(25);   file << "massU_broken[GeV]";
        file.width(25);   file << "massD_broken[GeV]";
        file.width(25);   file << "massS_broken[GeV]";
        file.width(30);   file << "energyDens_broken[GeVfm-3]";
        file.width(30);   file << "entropyDens_broken[fm-3]";
        file.width(25);   file << "sigmaU_broken[GeV]";
        file.width(25);   file << "sigmaD_broken[GeV]";
        file.width(25);   file << "sigmaS_broken[GeV]";
        file.width(25);   file << "rhoU_broken[fm-3]";
        file.width(25);   file << "rhoD_broken[fm-3]";
        file.width(25);   file << "rhoS_broken[fm-3]";

        file.width(25);   file << "rhoB_restored[fm-3]";
        file.width(25);   file << "massU_restored[GeV]";
        file.width(25);   file << "massD_restored[GeV]";
        file.width(25);   file << "massS_restored[GeV]";
        file.width(30);   file << "energyDens_restored[GeVfm-3]";
        file.width(30);   file << "entropyDens_restored[fm-3]";
        file.width(25);   file << "sigmaU_restored[GeV]";
        file.width(25);   file << "sigmaD_restored[GeV]";
        file.width(25);   file << "sigmaS_restored[GeV]";
        file.width(25);   file << "rhoU_restored[fm-3]";
        file.width(25);   file << "rhoD_restored[fm-3]";
        file.width(25);   file << "rhoS_restored[fm-3]";
        
        file << std::endl;
    }

    for (int i = 0; i < int(firstOrderLine.size()-1); i++)
    {
        ChiralTransitionPoint point = firstOrderLine[i];

        double T = point.temperature;
        SU3NJL3DCutoffParameters parametersNJL = point.parametersNJL;

        double mU_broken = point.upQuarkEffectiveMassBroken; 
        double mD_broken = point.downQuarkEffectiveMassBroken;
        double mS_broken = point.strangeQuarkEffectiveMassBroken;
        double effCP_broken = point.quarkEffectiveChemicalPotentialBroken;

        double mU_restored = point.upQuarkEffectiveMassRestored; 
        double mD_restored = point.downQuarkEffectiveMassRestored; 
        double mS_restored = point.strangeQuarkEffectiveMassRestored;
        double effCP_restored = point.quarkEffectiveChemicalPotentialRestored;

        //Broken phase - Calculate sigma and density for each quark flavour
        SU3NJL3DCutoffFixedTempRhoBEqualChemPot brokenPhase(
            parametersNJL,
            T,
            0.0,
            mU_broken,
            mD_broken,
            mS_broken,
            effCP_broken
        );
        brokenPhase.setSigmasDensitiesChemicalPotentials();
        double rhoU_broken = brokenPhase.getUpQuarkDensity();
        double rhoD_broken = brokenPhase.getDownQuarkDensity();
        double rhoS_broken = brokenPhase.getStrangeQuarkDensity();
        brokenPhase.setBaryonDensity( SU3BaryonDensity(rhoU_broken, rhoD_broken, rhoS_broken) );
        brokenPhase.setPressure( pressureVacuum );
        brokenPhase.setEnergyDensity( -pressureVacuum );
        brokenPhase.setEntropyDensity();

        //Restored phase - Calculate sigma and density for each quark flavour
        SU3NJL3DCutoffFixedTempRhoBEqualChemPot restoredPhase(
            parametersNJL,
            T,
            0.0,
            mU_restored,
            mD_restored,
            mS_restored,
            effCP_restored
        );
        restoredPhase.setSigmasDensitiesChemicalPotentials();
        double rhoU_restored = restoredPhase.getUpQuarkDensity();
        double rhoD_restored = restoredPhase.getDownQuarkDensity();
        double rhoS_restored = restoredPhase.getStrangeQuarkDensity();
        restoredPhase.setBaryonDensity( SU3BaryonDensity(rhoU_restored, rhoD_restored, rhoS_restored) );
        restoredPhase.setPressure( pressureVacuum );
        restoredPhase.setEnergyDensity( -pressureVacuum );
        restoredPhase.setEntropyDensity();

        // Print to file
        file.width(25);   file << effCP_broken;
        file.width(25);   file << brokenPhase.getBaryonChemicalPotential();
        file.width(25);   file << T;
        file.width(25);   file << brokenPhase.getPressure()/pow(hc_GeVfm,3);

        file.width(25);   file << brokenPhase.getBaryonDensity()/pow(hc_GeVfm,3);
        file.width(25);   file << brokenPhase.getUpQuarkEffectiveMass();
        file.width(25);   file << brokenPhase.getDownQuarkEffectiveMass();
        file.width(25);   file << brokenPhase.getStrangeQuarkEffectiveMass();
        file.width(30);   file << brokenPhase.getEnergyDensity()/pow(hc_GeVfm,3);
        file.width(30);   file << brokenPhase.getEntropyDensity()/pow(hc_GeVfm,3);
        file.width(25);   file << brokenPhase.getUpQuarkSigma();
        file.width(25);   file << brokenPhase.getDownQuarkSigma();
        file.width(25);   file << brokenPhase.getStrangeQuarkSigma();
        file.width(25);   file << brokenPhase.getUpQuarkDensity();
        file.width(25);   file << brokenPhase.getDownQuarkDensity();
        file.width(25);   file << brokenPhase.getStrangeQuarkDensity();

        file.width(25);   file << restoredPhase.getBaryonDensity()/pow(hc_GeVfm,3);
        file.width(25);   file << restoredPhase.getUpQuarkEffectiveMass();
        file.width(25);   file << restoredPhase.getDownQuarkEffectiveMass();
        file.width(25);   file << restoredPhase.getStrangeQuarkEffectiveMass();
        file.width(30);   file << restoredPhase.getEnergyDensity()/pow(hc_GeVfm,3);
        file.width(30);   file << restoredPhase.getEntropyDensity()/pow(hc_GeVfm,3);
        file.width(25);   file << restoredPhase.getUpQuarkSigma();
        file.width(25);   file << restoredPhase.getDownQuarkSigma();
        file.width(25);   file << restoredPhase.getStrangeQuarkSigma();
        file.width(25);   file << restoredPhase.getUpQuarkDensity();
        file.width(25);   file << restoredPhase.getDownQuarkDensity();
        file.width(25);   file << restoredPhase.getStrangeQuarkDensity();

        file << std::endl;
    }
}

void SU3NJL3DCutoffFixedTempRhoBEqualChemPot::evaluateFirstOrderLine(
    SU3NJL3DCutoffParameters& parameters,                                    
    double precisionVacuum,                                    
    MultiRootFindingMethod methodVacuum,                                    
    double upQuarkMassGuess, 
    double downQuarkMassGuess, 
    double strangeQuarkMassGuess,
    double minimumBaryonDensity_fmMinus3, 
    double maximumBaryonDensity_fmMinus3, 
    int numberOfPoints,
    double precisionZeroTempSol,
    MultiRootFindingMethod methodZeroTempSol,
    bool storeZeroTempSolToFile,
    double precisionTransitionPointSol, 
    MultiRootFindingMethod methodTransitionPointSol,
    double deltaT,
    double massDifferenceCEP
)
{
    SU3NJL3DCutoffVacuum vacuum = SU3NJL3DCutoffVacuum::evaluateVacuumMasses(
        parameters,                                    
        precisionVacuum,                                    
        methodVacuum,                                    
        upQuarkMassGuess, 
        downQuarkMassGuess, 
        strangeQuarkMassGuess
    );

    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint> firtOrderLine = 
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateFirstOrderLine(
        vacuum, 
        minimumBaryonDensity_fmMinus3*pow(hc_GeVfm,3), 
        maximumBaryonDensity_fmMinus3*pow(hc_GeVfm,3), 
        numberOfPoints, 
        precisionZeroTempSol, 
        methodZeroTempSol, 
        storeZeroTempSolToFile,
        precisionTransitionPointSol, 
        methodTransitionPointSol, 
        deltaT, 
        massDifferenceCEP
    );

    SU3NJL3DCutoffFixedTempRhoBEqualChemPot::writeToFile(
        vacuum, 
        firtOrderLine, 
        "SU3NJL3DCutoffFirstOrderLine_" + parameters.getParameterSetName() + ".dat",
        true
    );
}

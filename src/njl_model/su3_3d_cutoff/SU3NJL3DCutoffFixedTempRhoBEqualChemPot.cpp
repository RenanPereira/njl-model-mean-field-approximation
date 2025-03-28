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

SU3NJL3DCutoffFixedTempRhoBEqualChemPot::SU3NJL3DCutoffFixedTempRhoBEqualChemPot(SU3NJL3DCutoffParameters parametersNJLAux,
                                                                                 double temperatureAux)
{
	parametersNJL = parametersNJLAux;
    temperature = temperatureAux;
}

SU3NJL3DCutoffFixedTempRhoBEqualChemPot::SU3NJL3DCutoffFixedTempRhoBEqualChemPot(SU3NJL3DCutoffParameters parametersNJLAux,
                                                                                 double temperatureAux,
                                                                                 double baryonDensityAux,
                                                                                 double upQuarkEffectiveMassAux,
                                                                                 double downQuarkEffectiveMassAux,
                                                                                 double strangeQuarkEffectiveMassAux,
                                                                                 double quarkEffectiveChemicalPotentialAux)
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

    setTemperature(0.0); 
    setBaryonDensity(0.0);

    setUpQuarkEffectiveMass( vacuum.getUpQuarkEffectiveMass() );
	setDownQuarkEffectiveMass( vacuum.getDownQuarkEffectiveMass() );
	setStrangeQuarkEffectiveMass( vacuum.getStrangeQuarkEffectiveMass() );
    setQuarkEffectiveChemicalPotential(0.0);
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
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot solution(auxiliar);

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


bool SU3NJL3DCutoffFixedTempRhoBEqualChemPot::testSolution(double precision)
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


double SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculatePressure(double vacuumPressure)
{
    double pressureNJL = SU3NJL3DCutoffPressure(parametersNJL, temperature,
                                                upQuarkEffectiveMass, 
                                                downQuarkEffectiveMass, 
                                                strangeQuarkEffectiveMass, 
                                                quarkEffectiveChemicalPotential, 
                                                quarkEffectiveChemicalPotential, 
                                                quarkEffectiveChemicalPotential);

    pressureNJL = pressureNJL - vacuumPressure;

    return pressureNJL;
}


double SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateEnergyDensity(double vacuumEnergyDensity)
{
    double energyNJL = SU3NJL3DCutoffEnergyDensity(parametersNJL, temperature,
                                                   upQuarkEffectiveMass, 
                                                   downQuarkEffectiveMass, 
                                                   strangeQuarkEffectiveMass, 
                                                   quarkEffectiveChemicalPotential, 
                                                   quarkEffectiveChemicalPotential, 
                                                   quarkEffectiveChemicalPotential);

    energyNJL = energyNJL - vacuumEnergyDensity;

    return energyNJL;
}


double SU3NJL3DCutoffFixedTempRhoBEqualChemPot::calculateEntropyDensity()
{
    double entropyNJL = SU3NJL3DCutoffEntropyDensity(parametersNJL, temperature,
                                                     upQuarkEffectiveMass, 
                                                     downQuarkEffectiveMass, 
                                                     strangeQuarkEffectiveMass, 
                                                     quarkEffectiveChemicalPotential, 
                                                     quarkEffectiveChemicalPotential, 
                                                     quarkEffectiveChemicalPotential);

    return entropyNJL;
}


void SU3NJL3DCutoffFixedTempRhoBEqualChemPot::setSigmasAndDensitiesAndChemicalPotentials()
{
    //calculate the outputs of this solution
    NJLDimensionfulCouplings couplings = parametersNJL.getDimensionfulCouplings();

    NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();

    double Nc = parametersNJL.getNumberOfColours();
    
    double sigmaIntegralPrecision = parametersNJL.getSigmaIntegralPrecision();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

    double T = getTemperature();
    
    double mU = getUpQuarkEffectiveMass();
    double mD = getDownQuarkEffectiveMass();
    double mS = getStrangeQuarkEffectiveMass();
    double effCP = getQuarkEffectiveChemicalPotential();
   
    //Calculate the sigma field and density for each quark flavour
    double sigmaU = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCP, mU, sigmaIntegralPrecision);
    double sigmaD = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCP, mD, sigmaIntegralPrecision);
    double sigmaS = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effCP, mS, sigmaIntegralPrecision);

    double rhoU = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCP, mU, thermoIntegralPrecision);
    double rhoD = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCP, mD, thermoIntegralPrecision);
    double rhoS = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effCP, mS, thermoIntegralPrecision);

    //calculate quark chemical potential from the effective ones
    double cPU = SU3NJLQuarkChemicalPotential(couplings, effCP, rhoU, rhoD, rhoS, sigmaU, sigmaD, sigmaS);
    double cPD = SU3NJLQuarkChemicalPotential(couplings, effCP, rhoD, rhoS, rhoU, sigmaD, sigmaS, sigmaU);
    double cPS = SU3NJLQuarkChemicalPotential(couplings, effCP, rhoS, rhoU, rhoD, sigmaS, sigmaU, sigmaD);

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


vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> 
solveFromVacuumToFiniteBaryonDensity(SU3NJL3DCutoffVacuum vacuum, 
                                     double minimumBaryonDensity, double maximumBaryonDensity, int numberOfPoints, 
                                     double gapPrecision, MultiRootFindingMethod method)
{   
    // Analyse vacuum solution
    double pressureVacuum = vacuum.calculatePressure();
    double energyVacuum = vacuum.calculateEnergyDensity();
    cout << "pressureVacuum=" << pressureVacuum << "\n";
    cout << "energyDensityVacuum=" << energyVacuum << "\n";
    
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot vacuumSol(vacuum);
    vacuumSol.setPressure(0.0);
    vacuumSol.setEnergyDensity(0.0);
    vacuumSol.setEntropyDensity(0.0);

    // Build guesses from the vacuum solution
    double MuGuess = vacuumSol.getUpQuarkEffectiveMass();
    double MdGuess = vacuumSol.getDownQuarkEffectiveMass();
    double MsGuess = vacuumSol.getStrangeQuarkEffectiveMass();
    double effectiveCPGuess = vacuumSol.getUpQuarkEffectiveMass()*(1+1E-4);

    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> solutions;
    solutions.push_back(vacuumSol);
    double delta = (maximumBaryonDensity-minimumBaryonDensity)/(numberOfPoints-1);
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        SU3NJL3DCutoffFixedTempRhoBEqualChemPot inMediumSol(vacuum.getParametersNJL());
        inMediumSol.setTemperature(0.0); 
        inMediumSol.setBaryonDensity(0.0);

        //step in density
        double rhoB = minimumBaryonDensity + i*delta;
        inMediumSol.setBaryonDensity(rhoB);

        //find quark masses and effective chemical potential
        inMediumSol.solve(gapPrecision, method, MuGuess, MdGuess, MsGuess, effectiveCPGuess);

        //guesses for next step
        MuGuess = inMediumSol.getUpQuarkEffectiveMass();
        MdGuess = inMediumSol.getDownQuarkEffectiveMass();
        MsGuess = inMediumSol.getStrangeQuarkEffectiveMass();
        effectiveCPGuess = inMediumSol.getQuarkEffectiveChemicalPotential();  

        //print to console
        printf("rhoB=%.4f, Mu=%.4f, Md=%.4f, Ms=%.4f [GeV] \n", 
                rhoB/pow(hc_GeVfm, 3),
                inMediumSol.getUpQuarkEffectiveMass(), 
                inMediumSol.getDownQuarkEffectiveMass(), 
                inMediumSol.getStrangeQuarkEffectiveMass());
     

        //calculate thermodynamics
        double pressureMed = inMediumSol.calculatePressure(pressureVacuum);
        double energyMed = inMediumSol.calculateEnergyDensity(energyVacuum);
        double entropyMed = inMediumSol.calculateEntropyDensity();
        inMediumSol.setPressure( pressureMed );
        inMediumSol.setEnergyDensity( energyMed );
        inMediumSol.setEntropyDensity( entropyMed );

        //push to solutions vector if solution passes the test
        if (inMediumSol.testSolution(gapPrecision)==true){ solutions.push_back(inMediumSol); }
    }
    
    // Store calculation to file
    string filename = string("SU3NJL3DCutoffEqualChemPot_") + vacuum.getParametersNJL().getParameterSetName()
            	                                            + "T0.0"
                                                            + "rhoBMin" + trim0ToDot0(minimumBaryonDensity/pow(hc_GeVfm,3))
                                                            + "rhoBMax" + trim0ToDot0(maximumBaryonDensity/pow(hc_GeVfm,3))
                                                            + "N" + to_string(numberOfPoints)
                                                            + ".dat";
    writeSolutionsToFile(solutions, filename, true);

    return solutions;
}


void writeSolutionsToFile(vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> solutions, string fileName, bool columnsDescription)
{
    // Create file to store solution
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
        aux_file_NJL.width(25);   aux_file_NJL << 3.0*solutions[i].getQuarkEffectiveChemicalPotential();
   		aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getPressure()/pow(hc_GeVfm,3);
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getEnergyDensity()/pow(hc_GeVfm,3);
        aux_file_NJL.width(25);   aux_file_NJL << solutions[i].getEntropyDensity()/pow(hc_GeVfm,3);
		aux_file_NJL << std::endl;
    }
}


int SU3NJL3DCutoffEqualChemPotFixedTempChiralTransitionPoint(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{
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
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot brokenPhase(source.getParametersNJL(),
                                                        T,
                                                        0.0,
                                                        mU_broken,
                                                        mD_broken,
                                                        mS_broken,
                                                        effCP_broken);
    brokenPhase.setSigmasAndDensitiesAndChemicalPotentials();

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
    SU3NJL3DCutoffFixedTempRhoBEqualChemPot restoredPhase(source.getParametersNJL(),
                                                          T,
                                                          0.0,
                                                          mU_restored,
                                                          mD_restored,
                                                          mS_restored,
                                                          effCP_restored);
    restoredPhase.setSigmasAndDensitiesAndChemicalPotentials();

    double sigmaU_restored = restoredPhase.getUpQuarkSigma();
    double sigmaD_restored = restoredPhase.getDownQuarkSigma();
    double sigmaS_restored = restoredPhase.getStrangeQuarkSigma();

    //get all the necessary quantities
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


vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> findChiralTransitionPointsFixedTemperature(vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> solutions, double precision, MultiRootFindingMethod method)
{
    //run check that solutions have the same temperature
    for (int i = 0; i < int(solutions.size()); ++i)
    {   
        if ( fabs(solutions[1].getTemperature()-solutions[i].getTemperature())>0 ){ cout << "Temperature is not fixed in this vector of solutions!\n"; abort(); }
    }

    //finding guesses for transition point
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

    vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> transitionPoints;
    if ( first>0 && second>0 )
    {
        //estimate guesses
        double mU_broken = solutions[first].getUpQuarkEffectiveMass();
        double mD_broken = solutions[first].getDownQuarkEffectiveMass();
        double mS_broken = solutions[first].getStrangeQuarkEffectiveMass();
        double effCP_broken = solutions[first].getQuarkEffectiveChemicalPotential();
        double mU_restored = solutions[second].getUpQuarkEffectiveMass();
        double mD_restored = solutions[second].getDownQuarkEffectiveMass();
        double mS_restored = solutions[second].getStrangeQuarkEffectiveMass();
        double effCP_restored = solutions[second].getQuarkEffectiveChemicalPotential();

        //solve system
        double x[8];
        x[0] = mU_broken; 
        x[1] = mD_broken;  
        x[2] = mS_broken;  
        x[3] = 0.5*effCP_broken + 0.5*effCP_restored;
        x[4] = mU_restored; 
        x[5] = mD_restored;  
        x[6] = mS_restored;  
        x[7] = 0.5*effCP_broken + 0.5*effCP_restored;

        SU3NJL3DCutoffFixedTempRhoBEqualChemPot aux(solutions[0].getParametersNJL(), 
                                                    solutions[0].getTemperature());
        multiDimensionalRootFind(8, precision, &x[0], &aux, &SU3NJL3DCutoffEqualChemPotFixedTempChiralTransitionPoint, method);

        mU_broken = x[0]; 
        mD_broken = x[1];
        mS_broken = x[2]; 
        effCP_broken = x[3];
        mU_restored = x[4]; 
        mD_restored = x[5];  
        mS_restored = x[6];  
        effCP_restored = x[7]; 

        cout << "Broken phase:\n";
        cout << mU_broken << "\t" << mD_broken << "\t" << mS_broken << "\n";
        cout << effCP_broken << "\n";
        cout << "Restored phase:\n";
        cout << mU_restored << "\t" << mD_restored << "\t" << mS_restored << "\n";
        cout << effCP_restored << "\n";

    }

    return transitionPoints;
}
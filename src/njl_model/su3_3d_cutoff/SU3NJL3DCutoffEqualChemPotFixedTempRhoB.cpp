#include <cmath>
#include <iostream>
#include "physics_utils/physical_constants.h"
#include "njl_model/n_fermion_line_integrals/one_fermion_line_integral_3d_cutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffEqualChemPotFixedTempRhoB.h"
#include "utils/format_utils.h"

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


SU3NJL3DCutoffEqualChemPotFixedTempRhoB::SU3NJL3DCutoffEqualChemPotFixedTempRhoB(SU3NJL3DCutoffVacuum vacuum)
{   
    parametersNJL = vacuum.getParametersNJL();

    setTemperature(0.0); 
    setBaryonDensity(0.0);

    setUpQuarkEffectiveMass( vacuum.getUpQuarkEffectiveMass() );
	setDownQuarkEffectiveMass( vacuum.getDownQuarkEffectiveMass() );
	setStrangeQuarkEffectiveMass( vacuum.getStrangeQuarkEffectiveMass() );
    setQuarkEffectiveChemicalPotential(0.0);
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


double SU3NJL3DCutoffEqualChemPotFixedTempRhoB::calculateEnergyDensity(double vacuumEnergyDensity)
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


double SU3NJL3DCutoffEqualChemPotFixedTempRhoB::calculateEntropyDensity()
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


vector<SU3NJL3DCutoffEqualChemPotFixedTempRhoB> 
solveFromVacuumToFiniteBaryonDensity(SU3NJL3DCutoffVacuum vacuum, 
                                     double minimumBaryonDensity, double maximumBaryonDensity, int numberOfPoints, 
                                     double gapPrecision, MultiRootFindingMethod method)
{   
    // Analyse vacuum solution
    double pressureVacuum = vacuum.calculatePressure();
    double energyVacuum = vacuum.calculateEnergyDensity();
    cout << "pressureVacuum=" << pressureVacuum << "\n";
    cout << "energyDensityVacuum=" << energyVacuum << "\n";
    
    SU3NJL3DCutoffEqualChemPotFixedTempRhoB vacuumSol(vacuum);
    vacuumSol.setPressure(0.0);
    vacuumSol.setEnergyDensity(0.0);
    vacuumSol.setEntropyDensity(0.0);

    // Build guesses from the vacuum solution
    double MuGuess = vacuumSol.getUpQuarkEffectiveMass();
    double MdGuess = vacuumSol.getDownQuarkEffectiveMass();
    double MsGuess = vacuumSol.getStrangeQuarkEffectiveMass();
    double effectiveCPGuess = vacuumSol.getUpQuarkEffectiveMass()*(1+1E-4);

    vector<SU3NJL3DCutoffEqualChemPotFixedTempRhoB> solutions;
    solutions.push_back(vacuumSol);
    double delta = (maximumBaryonDensity-minimumBaryonDensity)/(numberOfPoints-1);
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        SU3NJL3DCutoffEqualChemPotFixedTempRhoB inMediumSol(vacuum.getParametersNJL());
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
        cout << "testSolution=" << inMediumSol.testSolution(gapPrecision) << "\n";
        cout << rhoB/pow(hc_GeVfm,3) << "\t" 
             << inMediumSol.getUpQuarkEffectiveMass() << "\t"
             << inMediumSol.getDownQuarkEffectiveMass() << "\t"
             << inMediumSol.getStrangeQuarkEffectiveMass() << "\n";  

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


void writeSolutionsToFile(vector<SU3NJL3DCutoffEqualChemPotFixedTempRhoB> solutions, string fileName, bool columnsDescription)
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

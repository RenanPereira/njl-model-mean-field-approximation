//Renan Câmara Pereira 2022-20??

#include <omp.h>
#include "command_line_processor.h"

using namespace std;


#include <gsl/gsl_complex_math.h>
#include "njl_model/njl_regularization_schemes.h"
#include "njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff.h"

std::string trim0ToDot0(const double );
void evaluateB0VSK0ToFile(const int , 
                          const double , 
                          const double ,
                          const NJL3DCutoffRegularizationScheme , 
                          const double ,
                          const double ,
                          const double ,
                          const double ,
                          const double ,
                          const double ,
                          const double ,
                          const double );

string trim0ToDot0(const double value) 
{
    string result = to_string(value);

    // Remove unnecessary trailing zeroes
    result.erase(result.find_last_not_of('0') + 1);

    // Ensure at least one digit remains after the decimal point (e.g., 0. -> 0.0)
    if (result.back() == '.') 
    {
        result += '0';
    }

    return result;
}

void evaluateB0VSK0ToFile(const int Npoints, 
                          const double k0LambdaRatioMin, 
                          const double k0LambdaRatioMax,
                          const NJL3DCutoffRegularizationScheme reguScheme, 
                          const double temperature,
                          const double effectiveChemicalPotential1,
                          const double effectiveChemicalPotential2,
                          const double threeMomentumCutoff,
                          const double effectiveMass1,
                          const double effectiveMass2,
                          const double threeMomentum,
                          const double integralPrecision)
{
    if (Npoints < 2 || k0LambdaRatioMax <= k0LambdaRatioMin) 
    {
        cout << "Error in evaluateB0VSK0ToFile: Npoints must be > 1 and k0LambdaRatioMax > k0LambdaRatioMin. Aborting.\n";
        abort();
    }

    vector<double> ratioExternalZeroMomentumTo3DCutoff(Npoints);
    vector<double> ReB0(Npoints);
    vector<double> ImB0(Npoints);
	double delta = ( k0LambdaRatioMax - k0LambdaRatioMin )/( Npoints - 1 );
	for (int i = 0; i < Npoints; ++i)
	{
		double k0 = ( k0LambdaRatioMin + i*delta )*threeMomentumCutoff;
        gsl_complex complexB0 = klevanskyB0Integral3DCutoff(reguScheme, 
                                                            temperature, 
                                                            effectiveChemicalPotential1, 
                                                            effectiveChemicalPotential2, 
                                                            threeMomentumCutoff, 
                                                            effectiveMass1, 
                                                            effectiveMass2, 
                                                            k0, 
                                                            threeMomentum, 
                                                            integralPrecision);

        ratioExternalZeroMomentumTo3DCutoff[i] = k0/threeMomentumCutoff;
        ReB0[i] = GSL_REAL(complexB0);
        ImB0[i] = GSL_IMAG(complexB0);
	}


    //Create file
    string filename = string("B0_vs_k0_") 
                    + "T"   + trim0ToDot0(temperature) 
    				+ "Cpi" + trim0ToDot0(effectiveChemicalPotential1) 
				    + "Cpj" + trim0ToDot0(effectiveChemicalPotential2)
    			    + "L"   + trim0ToDot0(threeMomentumCutoff)
    			    + "Mi"  + trim0ToDot0(effectiveMass1) 
    			    + "Mj"  + trim0ToDot0(effectiveMass2) 
    			    + "k"   + trim0ToDot0(threeMomentum) 
    			    + ".dat";
    ofstream fileB0;
    fileB0.open(filename, ofstream::out | ios::trunc);

    // Check if file is open successfully
    if (!fileB0.is_open()) 
    {
        cout << "Error: Unable to open file " << filename << endl;
        return;
    }

    fileB0.precision(15);
    fileB0.width(25);   fileB0 << "k03DCutoffRatio"; 
    fileB0.width(25);   fileB0 << "ReB0"; 
    fileB0.width(25);   fileB0 << "ImB0"; 
    fileB0 << endl;
    
    for (int i = 0; i < Npoints; ++i)
    {
        fileB0.width(25);   fileB0 << ratioExternalZeroMomentumTo3DCutoff[i]; 
        fileB0.width(25);   fileB0 << ReB0[i];
        fileB0.width(25);   fileB0 << ImB0[i]; 
        fileB0 << endl;
    }

    // Close the file explicitly
    fileB0.close();
}


int main(int argc, char* argv[])
{
	//START COUNTING TIME
    double start_s = omp_get_wtime();


    // Call the function and check its return value
	int processorResult = commandLineArgsProcessor(argc, argv);
    if ( processorResult==1 ) 
	{
        // If it returns 1, exit with failure code
        return 1;
    }
	else
	{
    	// Continue the main program execution if the function returns 0
    	std::cout << "\nCommands processed successfully, continuing execution..." << std::endl;
	}

/*
    NJL3DCutoffRegularizationScheme reguScheme = stringToNJL3DCutoffRegularizationScheme("CUTOFF_EVERYWHERE");
    double temperature = 0.0;
    double effectiveChemicalPotential1 = 0.0;
    double effectiveChemicalPotential2 = 0.0;
    double threeMomentumCutoff = 1.0;
    double effectiveMass1 = 0.4;
    double effectiveMass2 = 0.4;
    double threeMomentum = 0.0;
    double integralPrecision = 1E-8;
    evaluateB0VSK0ToFile(3000, 
                         -2.5, 
                         2.5,
                         reguScheme, 
                         temperature,
                         effectiveChemicalPotential1,
                         effectiveChemicalPotential2,
                         threeMomentumCutoff,
                         effectiveMass1,
                         effectiveMass2,
                         threeMomentum,
                         integralPrecision);
    
    threeMomentum = 0.5;
    evaluateB0VSK0ToFile(3000, 
                         -2.5, 
                         2.5,
                         reguScheme, 
                         temperature,
                         effectiveChemicalPotential1,
                         effectiveChemicalPotential2,
                         threeMomentumCutoff,
                         effectiveMass1,
                         effectiveMass2,
                         threeMomentum,
                         integralPrecision);

    threeMomentum = 1.0;
    evaluateB0VSK0ToFile(3000, 
                         -2.5, 
                         2.5,
                         reguScheme, 
                         temperature,
                         effectiveChemicalPotential1,
                         effectiveChemicalPotential2,
                         threeMomentumCutoff,
                         effectiveMass1,
                         effectiveMass2,
                         threeMomentum,
                         integralPrecision);

    threeMomentum = 1.5;
    evaluateB0VSK0ToFile(3000, 
                         -2.5, 
                         2.5,
                         reguScheme, 
                         temperature,
                         effectiveChemicalPotential1,
                         effectiveChemicalPotential2,
                         threeMomentumCutoff,
                         effectiveMass1,
                         effectiveMass2,
                         threeMomentum,
                         integralPrecision);
*/

/*
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
*/
/*
    //parameter set B
    double cutoff = 0.586967971572559;
    double gs = 3.0257845537123/pow(cutoff, 2);
    double kappa = -11.7520217701553/pow(cutoff, 5);
    double g1 = 36.5091272818204/pow(cutoff, 8);
    double g2 = -24.3337469126998/pow(cutoff, 8);
    double m0u = 0.00600659282357854;
    double m0d = 0.00600659282357854;
    double m0s = 0.138868437136382;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_SP8Q, gs, kappa, g1, g2);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setB");
*/

/*
    //parameter set C
    double cutoff = 0.586967971572559;
    double gs = 2.7596253718366/pow(cutoff, 2);
    double kappa = -11.7520217701553/pow(cutoff, 5);
    double g1 = 44.861215303826/pow(cutoff, 8);
    double g2 = -24.3337469126998/pow(cutoff, 8);
    double m0u = 0.00600659282357854;
    double m0d = 0.00600659282357854;
    double m0s = 0.138868437136382;

    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_SP8Q, gs, kappa, g1, g2);

    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("setC");


    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, HYBRIDS, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(1E-8) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";
*/

/*
    double chemPot = 0.318434158842783;
    double minimumTemperature = 0.040;
    double maximumTemperature = 0.300;
    int numberOfPointsFromVacToMinTemp = 200; 
    int numberOfPointsMinTempToChemPot = 200;
    int numberOfPointsFromMinToMaxTemp = 261;
    bool largeAngleScatteringContribution = false;
    IntegratedCrossSectionApproximationMethod approximationMethod = completeCOV;
    double propagatorIntegralPrecision = 1E-7;
    double crossSectionIntegralPrecision = 1E-4;
    double integratedCrossSectionIntegralPrecision_dXdY = 1E-12;
    double integratedCrossSectionIntegralPrecision_dX = 1E-3;
    evaluateIntegratedCrossSectionsWithFixedChemicalPotential(vacuum,
                                                              1E-8,
                                                              chemPot,
                                                              minimumTemperature, 
                                                              maximumTemperature, 
                                                              numberOfPointsFromVacToMinTemp, 
                                                              numberOfPointsMinTempToChemPot,
                                                              numberOfPointsFromMinToMaxTemp, 
                                                              largeAngleScatteringContribution, 
                                                              approximationMethod,
                                                              propagatorIntegralPrecision,
                                                              crossSectionIntegralPrecision,
                                                              integratedCrossSectionIntegralPrecision_dXdY,
                                                              integratedCrossSectionIntegralPrecision_dX);

    largeAngleScatteringContribution = true;
    evaluateIntegratedCrossSectionsWithFixedChemicalPotential(vacuum,
                                                              1E-8,
                                                              chemPot,
                                                              minimumTemperature, 
                                                              maximumTemperature, 
                                                              numberOfPointsFromVacToMinTemp, 
                                                              numberOfPointsMinTempToChemPot,
                                                              numberOfPointsFromMinToMaxTemp, 
                                                              largeAngleScatteringContribution, 
                                                              approximationMethod,
                                                              propagatorIntegralPrecision,
                                                              crossSectionIntegralPrecision,
                                                              integratedCrossSectionIntegralPrecision_dXdY,
                                                              integratedCrossSectionIntegralPrecision_dX);
*/
/*
    double minimumTemperature = 0.120;
    double maximumTemperature = 0.300;
    int numberOfPointsFromVacToMinTemp = 200; 
    int numberOfPointsFromMinToMaxTemp = 181;
    bool largeAngleScatteringContribution = false;
    IntegratedCrossSectionApproximationMethod approximationMethod = completeCOV;
    double propagatorIntegralPrecision = 1E-7;
    double crossSectionIntegralPrecision = 1E-4;
    double integratedCrossSectionIntegralPrecision_dXdY = 1E-12;
    double integratedCrossSectionIntegralPrecision_dX = 1E-3;
    evaluateIntegratedCrossSectionsWithZeroChemicalPotential(vacuum,
                                                             1E-8,
                                                             minimumTemperature, 
                                                             maximumTemperature, 
                                                             numberOfPointsFromVacToMinTemp, 
                                                             numberOfPointsFromMinToMaxTemp, 
                                                             largeAngleScatteringContribution, 
                                                             approximationMethod,
                                                             propagatorIntegralPrecision,
                                                             crossSectionIntegralPrecision,
                                                             integratedCrossSectionIntegralPrecision_dXdY,
                                                             integratedCrossSectionIntegralPrecision_dX);

    largeAngleScatteringContribution = true;
    evaluateIntegratedCrossSectionsWithZeroChemicalPotential(vacuum,
                                                             1E-8,
                                                             minimumTemperature, 
                                                             maximumTemperature, 
                                                             numberOfPointsFromVacToMinTemp, 
                                                             numberOfPointsFromMinToMaxTemp, 
                                                             largeAngleScatteringContribution, 
                                                             approximationMethod,
                                                             propagatorIntegralPrecision,
                                                             crossSectionIntegralPrecision,
                                                             integratedCrossSectionIntegralPrecision_dXdY,
                                                             integratedCrossSectionIntegralPrecision_dX);
*/




/*
/////////////////////////////////////////////////////////////////////////////////////
//Neutron star equation of state stuff


    //parameter set (Renan Master thesis parameter set)
    double cutoff = 0.6023;
    double gs = 2*( 1.835/pow(cutoff,2) );
    double kappa = -12.360/pow(cutoff,5);
    double m0u = 0.0055;
    double m0d = 0.0055;
    double m0s = 0.1407;

    //double gOmega1 = 0.2*pow(0.5*gs, 1);
    //double gOmega2 = 5.0*pow(0.5*gs, 4);
    //double gOmega3 = -10.0*pow(0.5*gs, 7);
    //double gOmega4 = 7.0*pow(0.5*gs, 10);

    //Fix Lagrangian dimensionful couplings
    //NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_VP4Q_VP8Q_VP12Q_VP16Q, gs, kappa, gOmega1, gOmega2, gOmega3, gOmega4);

    vector<double> gOmegaAdimensional = {0.8, 1.0, -3.0, 3.0, -1.0};
    vector<double> gOmegaDimensionful = multiQuarkVPCouplingWithDimensions(gOmegaAdimensional, 0.5*gs);
    
    //Fix Lagrangian dimensionful couplings
    NJLDimensionfulCouplings couplings(SP4Q_DET2NFQ_VPMULTIQ, gs, kappa, gOmegaDimensionful);


    //Create NJL parameter set
    SU3NJL3DCutoffParameters parameters(CUTOFF_EVERYWHERE, cutoff, couplings, m0u, m0d, m0s);
    parameters.setParameterSetName("renanMasterThesis");


    //solve model in the vacuum
    double gapPrecision = 1E-8;
    SU3NJL3DCutoffVacuum vacuum(parameters);
    vacuum.solve(gapPrecision, HYBRIDS, 0.3, 0.3, 0.5);

    cout << "Vacuum effective masses: \n";
    cout << "testSolution=" << vacuum.testSolution(gapPrecision) << "\n";
    cout << "Mu=" << vacuum.getUpQuarkEffectiveMass() << "GeV" << "\t" 
         << "Md=" << vacuum.getDownQuarkEffectiveMass() << "GeV" << "\t" 
         << "Ms=" << vacuum.getStrangeQuarkEffectiveMass() << "GeV" << "\n";


    double rhoi = 1E-5*pow(hc_GeVfm, 3);
    double rhof = 3.80*pow(hc_GeVfm, 3);
    int NrhoB = 8000;
    writeBetaEquilibriumEOSAtZeroTemperatureToFile(vacuum, rhoi, rhof, NrhoB, gapPrecision, HYBRIDS, "eos.dat");


    std::ofstream file;
    file.open("gOmega.dat", std::fstream::in | std::ofstream::out | std::ios::trunc);
    string gOmegaId;
    for (int i = 0; i < int( gOmegaAdimensional.size() ); ++i)
    {
        gOmegaId = gOmegaId + "_" + to_string(gOmegaAdimensional[i]);
    }
    file << gOmegaId;
    file.close();
*/


    //STOP CLOCK AND PRINT RUN TIME
    double stop_s = omp_get_wtime();
    double run_time = (stop_s-start_s);
    std::cout << "Run Time: " << run_time << std::endl;

	return 0;
}




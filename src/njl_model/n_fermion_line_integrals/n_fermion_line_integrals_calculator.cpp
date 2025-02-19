#include "njl_model/n_fermion_line_integrals/n_fermion_line_integrals_calculator.h"
#include "njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff.h"
#include "utils/format_utils.h"
#include <fstream>

using namespace std;


void evaluateKlevanskyB0Integral3DCutoffVsZeroMomentumToFile(
    const int Npoints, 
    const double k0ToLambdaRatioMin, 
    const double k0ToLambdaRatioMax,
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
    if (Npoints < 2 || k0ToLambdaRatioMax <= k0ToLambdaRatioMin) 
    {
        cout << "Error in evaluateB0VSK0ToFile: Npoints must be > 1 and k0ToLambdaRatioMax > k0ToLambdaRatioMin. Aborting.\n";
        abort();
    }

    vector<double> ratioZeroMomentumTo3DCutoff(Npoints);
    vector<double> RealB0(Npoints);
    vector<double> ImagB0(Npoints);
	double delta = ( k0ToLambdaRatioMax - k0ToLambdaRatioMin )/( Npoints - 1 );
	for (int i = 0; i < Npoints; ++i)
	{
		double k0 = ( k0ToLambdaRatioMin + i*delta )*threeMomentumCutoff;
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

        ratioZeroMomentumTo3DCutoff[i] = k0/threeMomentumCutoff;
        RealB0[i] = GSL_REAL(complexB0);
        ImagB0[i] = GSL_IMAG(complexB0);
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
    fileB0.width(25);   fileB0 << "k0To3DCutoffRatio"; 
    fileB0.width(25);   fileB0 << "RealB0"; 
    fileB0.width(25);   fileB0 << "ImagB0"; 
    fileB0 << endl;
    
    for (int i = 0; i < Npoints; ++i)
    {
        fileB0.width(25);   fileB0 << ratioZeroMomentumTo3DCutoff[i]; 
        fileB0.width(25);   fileB0 << RealB0[i];
        fileB0.width(25);   fileB0 << ImagB0[i]; 
        fileB0 << endl;
    }

    // Close the file explicitly
    fileB0.close();
}

void evaluateKlevanskyB0Integral3DCutoffVsThreeMomentumToFile(
    const int Npoints, 
    const double kToLambdaRatioMin, 
    const double kToLambdaRatioMax,
    const NJL3DCutoffRegularizationScheme reguScheme, 
    const double temperature,
    const double effectiveChemicalPotential1,
    const double effectiveChemicalPotential2,
    const double threeMomentumCutoff,
    const double effectiveMass1,
    const double effectiveMass2,
    const double k0,
    const double integralPrecision)
{
    if (Npoints < 2 || kToLambdaRatioMax <= kToLambdaRatioMin) 
    {
        cout << "Error in evaluateKlevanskyB0IntegralVsThreeMomentumToFile: Npoints must be > 1 and kToLambdaRatioMax > k0ToLambdaRatioMin. Aborting.\n";
        abort();
    }

    vector<double> kToLambda(Npoints);
    vector<double> RealB0(Npoints);
    vector<double> ImagB0(Npoints);
	double delta = ( kToLambdaRatioMax - kToLambdaRatioMin )/( Npoints - 1 );
	for (int i = 0; i < Npoints; ++i)
	{
		double k = ( kToLambdaRatioMin + i*delta )*threeMomentumCutoff;
        gsl_complex complexB0 = klevanskyB0Integral3DCutoff(reguScheme, 
                                                            temperature, 
                                                            effectiveChemicalPotential1, 
                                                            effectiveChemicalPotential2, 
                                                            threeMomentumCutoff, 
                                                            effectiveMass1, 
                                                            effectiveMass2, 
                                                            k0, 
                                                            k, 
                                                            integralPrecision);

        kToLambda[i] = k/threeMomentumCutoff;
        RealB0[i] = GSL_REAL(complexB0);
        ImagB0[i] = GSL_IMAG(complexB0);
	}


    //Create file
    string filename = string("B0_vs_k_") 
                    + "T"   + trim0ToDot0(temperature) 
    				+ "Cpi" + trim0ToDot0(effectiveChemicalPotential1) 
				    + "Cpj" + trim0ToDot0(effectiveChemicalPotential2)
    			    + "L"   + trim0ToDot0(threeMomentumCutoff)
    			    + "Mi"  + trim0ToDot0(effectiveMass1) 
    			    + "Mj"  + trim0ToDot0(effectiveMass2) 
    			    + "k0"  + trim0ToDot0(k0) 
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
    fileB0.width(25);   fileB0 << "kTo3DCutoffRatio"; 
    fileB0.width(25);   fileB0 << "RealB0"; 
    fileB0.width(25);   fileB0 << "ImagB0"; 
    fileB0 << endl;
    
    for (int i = 0; i < Npoints; ++i)
    {
        fileB0.width(25);   fileB0 << kToLambda[i]; 
        fileB0.width(25);   fileB0 << RealB0[i];
        fileB0.width(25);   fileB0 << ImagB0[i]; 
        fileB0 << endl;
    }

    // Close the file explicitly
    fileB0.close();
}

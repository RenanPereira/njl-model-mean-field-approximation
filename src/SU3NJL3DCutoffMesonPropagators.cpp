#include <cmath>
#include <iostream>
#include <gsl/gsl_complex_math.h>
#include "OneFermionLineIntegral.h"
#include "TwoFermionLineIntegral.h"
#include "SU3NJL3DCutoffMesonProjectors.h"
#include "SU3NJL3DCutoffMesonPropagators.h"

using namespace std;


//In the scalar and pseudoscalar polarization operators a finite width is included with the following recipe: k0 -> k0 - iGamma/2



gsl_complex pseudoscalarPolarizationOperator3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double Nc, 
													 double T, double effCP1, double effCP2, double M1, double M2, 
													 double k0, double k, double Gamma, double integralPrecision)
{	
	double deltaReal = pow((M1-M2), 2) - pow((k0 + effCP1 - effCP2), 2) + pow(k,2) + pow(Gamma,2)/4.0;
	double deltaImag = Gamma*(k0 + effCP1 - effCP2);
	gsl_complex delta = gsl_complex_rect(deltaReal, deltaImag);

	//calculate contributions coming from the one fermion line integrals
	gsl_complex klevA1 = klevanskyAIntegral3DCutoff(reguScheme, cutoff, T, effCP1, M1, k, integralPrecision);
    gsl_complex klevA2 = klevanskyAIntegral3DCutoff(reguScheme, cutoff, T, effCP2, M2, k, integralPrecision);

    //calculate contribution coming from the two fermion line integral
    gsl_complex klevB0 = klevanskyB0Integral3DCutoff(reguScheme, T, effCP1, effCP2, cutoff, M1, M2, k0, k, integralPrecision);

    //build the pseudoscalar quark polarization operator
    gsl_complex pseudoscalarPol;
    pseudoscalarPol = gsl_complex_add(klevA1, klevA2);
    pseudoscalarPol = gsl_complex_add(pseudoscalarPol, gsl_complex_mul(klevB0, delta));
    pseudoscalarPol = gsl_complex_mul_real(pseudoscalarPol, -2.0*Nc/( 16.0*(M_PI*M_PI) ) );

    return pseudoscalarPol;
}


gsl_complex pseudoscalarPolarizationOperator3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double Nc, 
												     double T, double effCP1, double effCP2, double M1, double M2, 
												     double k0, double k, double integralPrecision)
{
	gsl_complex pseudoscalarPol = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effCP1, effCP2, M1, M2, k0, k, 0.0, integralPrecision);

	return pseudoscalarPol;
}


gsl_complex scalarPolarizationOperator3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double Nc, 
											   double T, double effCP1, double effCP2, double M1, double M2, 
											   double k0, double k, double Gamma, double integralPrecision)
{
	double deltaReal = pow((M1+M2), 2) - pow((k0 + effCP1 - effCP2), 2) + pow(k,2) + pow(Gamma,2)/4.0;
    double deltaImag = Gamma*(k0 + effCP1 - effCP2);
    gsl_complex delta = gsl_complex_rect(deltaReal, deltaImag);

	//calculate contributions coming from the one fermion line integrals
	gsl_complex klevA1 = klevanskyAIntegral3DCutoff(reguScheme, cutoff, T, effCP1, M1, k, integralPrecision);
    gsl_complex klevA2 = klevanskyAIntegral3DCutoff(reguScheme, cutoff, T, effCP2, M2, k, integralPrecision);

    //calculate contribution coming from the two fermion line integral
    gsl_complex klevB0 = klevanskyB0Integral3DCutoff(reguScheme, T, effCP1, effCP2, cutoff, M1, M2, k0, k, integralPrecision);

	//build the pseudoscalar quark polarization operator
    gsl_complex scalarPol;
    scalarPol = gsl_complex_add(klevA1, klevA2);
    scalarPol = gsl_complex_add(scalarPol, gsl_complex_mul(klevB0, delta));
    scalarPol = gsl_complex_mul_real(scalarPol, -2.0*Nc/( 16.0*(M_PI*M_PI) ) );

    return scalarPol;
}


gsl_complex scalarPolarizationOperator3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double Nc, 
											   double T, double effCP1, double effCP2, double M1, double M2, 
											   double k0, double k, double integralPrecision)
{
	gsl_complex scalarPol = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effCP1, effCP2, M1, M2, k0, k, 0.0, integralPrecision);

	return scalarPol;
}


////////////////////////////////////////////////////////////////////////////////////////
//Pseudoscalar Meson Propagators


gsl_complex pionPlusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		   double effMassU, double effMassD, double effMassS, 
					  		   double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double P11 = projector.pseudoscalar11();

	//calculate appropriate polarization function
	gsl_complex pseudoscalarPolUD = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotU, effChemPotD, effMassU, effMassD, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(pseudoscalarPolUD, -4.0*P11);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*P11, 0.0);

    gsl_complex pionPropagator = gsl_complex_div(numerator, denominator);

    return pionPropagator;
}


gsl_complex pionPlusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		   double effMassU, double effMassD, double effMassS, 
					  		   double k0, double k, double integralPrecision)
{
	gsl_complex pionPropagator = 
	pionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return pionPropagator;
}


gsl_complex pionMinusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		    double effChemPotU, double effChemPotD, double effChemPotS, 
					  		    double effMassU, double effMassD, double effMassS, 
					  		    double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double P11 = projector.pseudoscalar11();

	//calculate appropriate polarization function
	gsl_complex pseudoscalarPolDU = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotD, effChemPotU, effMassD, effMassU, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(pseudoscalarPolDU, -4.0*P11);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*P11, 0.0);

    gsl_complex pionPropagator = gsl_complex_div(numerator, denominator);

    return pionPropagator;
}


gsl_complex pionMinusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		    double effChemPotU, double effChemPotD, double effChemPotS, 
					  		    double effMassU, double effMassD, double effMassS, 
					  		    double k0, double k, double integralPrecision)
{
	gsl_complex pionPropagator = 
	pionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return pionPropagator;
}


gsl_complex kaonPlusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		   double effMassU, double effMassD, double effMassS, 
					  		   double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double P44 = projector.pseudoscalar44();

	//calculate appropriate polarization function
	gsl_complex pseudoscalarPolUS = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotU, effChemPotS, effMassU, effMassS, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(pseudoscalarPolUS, -4.0*P44);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*P44, 0.0);

    gsl_complex kaonPropagator = gsl_complex_div(numerator, denominator);

    return kaonPropagator;
}


gsl_complex kaonPlusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		   double effMassU, double effMassD, double effMassS, 
					  		   double k0, double k, double integralPrecision)
{
	gsl_complex kaonPropagator = 
	kaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return kaonPropagator;
}


gsl_complex kaonMinusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		    double effChemPotU, double effChemPotD, double effChemPotS, 
					  		    double effMassU, double effMassD, double effMassS, 
					  		    double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double P44 = projector.pseudoscalar44();

	//calculate appropriate polarization function
	gsl_complex pseudoscalarPolSU = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotS, effChemPotU, effMassS, effMassU, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(pseudoscalarPolSU, -4.0*P44);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*P44, 0.0);

    gsl_complex kaonPropagator = gsl_complex_div(numerator, denominator);

    return kaonPropagator;
}


gsl_complex kaonMinusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		    double effChemPotU, double effChemPotD, double effChemPotS, 
					  		    double effMassU, double effMassD, double effMassS, 
					  		    double k0, double k, double integralPrecision)
{
	gsl_complex kaonPropagator = 
	kaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return kaonPropagator;
}


gsl_complex neutralKaonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		      double effChemPotU, double effChemPotD, double effChemPotS, 
					  		      double effMassU, double effMassD, double effMassS, 
					  		      double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double P66 = projector.pseudoscalar66();

	//calculate appropriate polarization function
	gsl_complex pseudoscalarPolDS = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotD, effChemPotS, effMassD, effMassS, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(pseudoscalarPolDS, -4.0*P66);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*P66, 0.0);

    gsl_complex kaonPropagator = gsl_complex_div(numerator, denominator);

    return kaonPropagator;
}


gsl_complex neutralKaonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		      double effChemPotU, double effChemPotD, double effChemPotS, 
					  		      double effMassU, double effMassD, double effMassS, 
					  		      double k0, double k, double integralPrecision)
{
	gsl_complex kaonPropagator = 
	neutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return kaonPropagator;
}


gsl_complex antiNeutralKaonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          double effMassU, double effMassD, double effMassS, 
					  		          double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double P66 = projector.pseudoscalar66();

	//calculate appropriate polarization function
	gsl_complex pseudoscalarPolSD = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotS, effChemPotD, effMassS, effMassD, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(pseudoscalarPolSD, -4.0*P66);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*P66, 0.0);

    gsl_complex kaonPropagator = gsl_complex_div(numerator, denominator);

    return kaonPropagator;
}


gsl_complex antiNeutralKaonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          double effMassU, double effMassD, double effMassS, 
					  		          double k0, double k, double integralPrecision)
{
	gsl_complex kaonPropagator = 
	antiNeutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return kaonPropagator;
}


ComplexSquareMatrixGSL neutral038PseudoscalarsPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          				     double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          				     double effMassU, double effMassD, double effMassS, 
					  		          			         double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//Set matrix with meson projection operators
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double P00 = projector.pseudoscalar00();
	double P03 = projector.pseudoscalar03();
	double P08 = projector.pseudoscalar08();
	double P33 = projector.pseudoscalar33();
	double P38 = projector.pseudoscalar38();
	double P88 = projector.pseudoscalar88();

	ComplexSquareMatrixGSL pAB(3);

	pAB.setValue(0, 0, gsl_complex_rect(2.0*P00, 0.0));
	pAB.setValue(0, 1, gsl_complex_rect(2.0*P03, 0.0));
	pAB.setValue(0, 2, gsl_complex_rect(2.0*P08, 0.0));

	pAB.setValue(1, 0, gsl_complex_rect(2.0*P03, 0.0));
	pAB.setValue(1, 1, gsl_complex_rect(2.0*P33, 0.0));
	pAB.setValue(1, 2, gsl_complex_rect(2.0*P38, 0.0));

	pAB.setValue(2, 0, gsl_complex_rect(2.0*P08, 0.0));
	pAB.setValue(2, 1, gsl_complex_rect(2.0*P38, 0.0));
	pAB.setValue(2, 2, gsl_complex_rect(2.0*P88, 0.0));

	//Set matrix with quark polarization operators
	gsl_complex pseudoscalarPolUU = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotU, effChemPotU, effMassU, effMassU, k0, k, Gamma, integralPrecision);
	gsl_complex pseudoscalarPolDD = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotD, effChemPotD, effMassD, effMassD, k0, k, Gamma, integralPrecision);
	gsl_complex pseudoscalarPolSS = 
	pseudoscalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotS, effChemPotS, effMassS, effMassS, k0, k, Gamma, integralPrecision);

	//PI00 = (2./3.)*( PIuu + PIdd + PIss)
	gsl_complex pseudoscalarPol00;
	pseudoscalarPol00 = gsl_complex_add(pseudoscalarPolUU, pseudoscalarPolDD);
    pseudoscalarPol00 = gsl_complex_add(pseudoscalarPol00, pseudoscalarPolSS);
    pseudoscalarPol00 = gsl_complex_mul_real(pseudoscalarPol00, 2.0/3.0);

    //PI03 = sqrt(2./3.)*( PIuu - PIdd )
    gsl_complex pseudoscalarPol03;
    pseudoscalarPol03 = gsl_complex_sub(pseudoscalarPolUU, pseudoscalarPolDD);
    pseudoscalarPol03 = gsl_complex_mul_real(pseudoscalarPol03, sqrt(2.0/3.0));

    //PI08 = (sqrt(2)/3.)*( PIuu + PIdd - 2 PIss )
    gsl_complex pseudoscalarPol08;
    pseudoscalarPol08 = gsl_complex_add(pseudoscalarPolUU, pseudoscalarPolDD);
    pseudoscalarPol08 = gsl_complex_add(pseudoscalarPol08, gsl_complex_mul_real(pseudoscalarPolSS, -2.0));
    pseudoscalarPol08 = gsl_complex_mul_real(pseudoscalarPol08, sqrt(2.0)/3.0);

	//PI33 = ( PIuu + PIdd )
    gsl_complex pseudoscalarPol33;
    pseudoscalarPol33 = gsl_complex_add(pseudoscalarPolUU, pseudoscalarPolDD); 

    //PI38 = (1/sqrt(3.))*( PIuu - PIdd )
    gsl_complex pseudoscalarPol38;
    pseudoscalarPol38 = gsl_complex_sub(pseudoscalarPolUU, pseudoscalarPolDD);
    pseudoscalarPol38 = gsl_complex_mul_real(pseudoscalarPol38, 1.0/sqrt(3.0));

    //PI88 = (1./3.)*( Piuu + PIdd + 4 PIss )
    gsl_complex pseudoscalarPol88;
    pseudoscalarPol88 = gsl_complex_add(pseudoscalarPolUU, pseudoscalarPolDD);
    pseudoscalarPol88 = gsl_complex_add(pseudoscalarPol88, gsl_complex_mul_real(pseudoscalarPolSS, 4.0));
    pseudoscalarPol88 = gsl_complex_mul_real(pseudoscalarPol88, 1.0/3.0);

    ComplexSquareMatrixGSL pseudoscalarPol(3);

    pseudoscalarPol.setValue(0, 0, pseudoscalarPol00);
    pseudoscalarPol.setValue(0, 1, pseudoscalarPol03);
    pseudoscalarPol.setValue(0, 2, pseudoscalarPol08);

    pseudoscalarPol.setValue(1, 0, pseudoscalarPol03);
    pseudoscalarPol.setValue(1, 1, pseudoscalarPol33);
    pseudoscalarPol.setValue(1, 2, pseudoscalarPol38);

    pseudoscalarPol.setValue(2, 0, pseudoscalarPol08);
    pseudoscalarPol.setValue(2, 1, pseudoscalarPol38);
    pseudoscalarPol.setValue(2, 2, pseudoscalarPol88);

    //Set inverse meson propagator
    ComplexSquareMatrixGSL inverseNeutral038Propagator = subtract(pAB.inverse(), pseudoscalarPol);

    //return propagator, invert the inverse propagator!
	return inverseNeutral038Propagator.inverse();
}


ComplexSquareMatrixGSL neutral038PseudoscalarsPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          				     double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          				     double effMassU, double effMassD, double effMassS, 
					  		          				     double k0, double k, double integralPrecision)
{
	ComplexSquareMatrixGSL neutral038Propagator = 
	neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return neutral038Propagator;
}


////////////////////////////////////////////////////////////////////////////////////////
//Scalar Meson Propagators


gsl_complex sigmaPionPlusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		   		double effChemPotU, double effChemPotD, double effChemPotS, 
					  		   		double effMassU, double effMassD, double effMassS, 
					  		   		double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double S11 = projector.scalar11();

	//calculate appropriate polarization function
	gsl_complex scalarPolUD = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotU, effChemPotD, effMassU, effMassD, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(scalarPolUD, -4.0*S11);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*S11, 0.0);

    gsl_complex sigmaPionPropagator = gsl_complex_div(numerator, denominator);

    return sigmaPionPropagator;
}


gsl_complex sigmaPionPlusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		   		double effChemPotU, double effChemPotD, double effChemPotS, 
					  		   		double effMassU, double effMassD, double effMassS, 
					  		   		double k0, double k, double integralPrecision)
{
	gsl_complex sigmaPionPropagator = 
	sigmaPionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return sigmaPionPropagator;
}


gsl_complex sigmaPionMinusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		    	 double effChemPotU, double effChemPotD, double effChemPotS, 
					  		    	 double effMassU, double effMassD, double effMassS, 
					  		    	 double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double S11 = projector.scalar11();

	//calculate appropriate polarization function
	gsl_complex scalarPolDU = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotD, effChemPotU, effMassD, effMassU, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(scalarPolDU, -4.0*S11);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*S11, 0.0);

    gsl_complex sigmaPionPropagator = gsl_complex_div(numerator, denominator);


    return sigmaPionPropagator;
}


gsl_complex sigmaPionMinusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		    	 double effChemPotU, double effChemPotD, double effChemPotS, 
					  		    	 double effMassU, double effMassD, double effMassS, 
					  		    	 double k0, double k, double integralPrecision)
{
	gsl_complex sigmaPionPropagator = 
	sigmaPionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return sigmaPionPropagator;
}


gsl_complex sigmaKaonPlusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		   		double effChemPotU, double effChemPotD, double effChemPotS, 
					  		   		double effMassU, double effMassD, double effMassS, 
					  		   		double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double S44 = projector.scalar44();

	//calculate appropriate polarization function
	gsl_complex scalarPolUS = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotU, effChemPotS, effMassU, effMassS, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(scalarPolUS, -4.0*S44);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*S44, 0.0);

    gsl_complex sigmaKaonPropagator = gsl_complex_div(numerator, denominator);

    return sigmaKaonPropagator;
}


gsl_complex sigmaKaonPlusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		   		double effChemPotU, double effChemPotD, double effChemPotS, 
					  		   		double effMassU, double effMassD, double effMassS, 
					  		   		double k0, double k, double integralPrecision)
{
	gsl_complex sigmaKaonPropagator = 
	sigmaKaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return sigmaKaonPropagator;
}


gsl_complex sigmaKaonMinusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		    	 double effChemPotU, double effChemPotD, double effChemPotS, 
					  		    	 double effMassU, double effMassD, double effMassS, 
					  		    	 double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double S44 = projector.scalar44();

	//calculate appropriate polarization function
	gsl_complex scalarPolSU = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotS, effChemPotU, effMassS, effMassU, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(scalarPolSU, -4.0*S44);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*S44, 0.0);

    gsl_complex sigmaKaonPropagator = gsl_complex_div(numerator, denominator);

    return sigmaKaonPropagator;
}


gsl_complex sigmaKaonMinusPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		    	 double effChemPotU, double effChemPotD, double effChemPotS, 
					  		    	 double effMassU, double effMassD, double effMassS, 
					  		    	 double k0, double k, double integralPrecision)
{
	gsl_complex sigmaKaonPropagator = 
	sigmaKaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return sigmaKaonPropagator;
}


gsl_complex neutralSigmaKaonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		      	   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		      	   double effMassU, double effMassD, double effMassS, 
					  		      	   double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double S66 = projector.scalar66();

	//calculate appropriate polarization function
	gsl_complex scalarPolDS = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotD, effChemPotS, effMassD, effMassS, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(scalarPolDS, -4.0*S66);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*S66, 0.0);

    gsl_complex sigmaKaonPropagator = gsl_complex_div(numerator, denominator);

    return sigmaKaonPropagator;
}


gsl_complex neutralSigmaKaonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		      	   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		      	   double effMassU, double effMassD, double effMassS, 
					  		      	   double k0, double k, double integralPrecision)
{
	gsl_complex sigmaKaonPropagator = 
	neutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return sigmaKaonPropagator;
}


gsl_complex antiNeutralSigmaKaonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          	   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          	   double effMassU, double effMassD, double effMassS, 
					  		          	   double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//calculate appropriate projector
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double S66 = projector.scalar66();

	//calculate appropriate polarization function
	gsl_complex scalarPolSD = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotS, effChemPotD, effMassS, effMassD, k0, k, Gamma, integralPrecision);

	//build propagator
	gsl_complex denominator, numerator;
    denominator = gsl_complex_mul_real(scalarPolSD, -4.0*S66);
    denominator = gsl_complex_add_real(denominator, 1.0);
    numerator = gsl_complex_rect(2.0*S66, 0.0);

    gsl_complex sigmaKaonPropagator = gsl_complex_div(numerator, denominator);

    return sigmaKaonPropagator;
}


gsl_complex antiNeutralSigmaKaonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          	   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          	   double effMassU, double effMassD, double effMassS, 
					  		          	   double k0, double k, double integralPrecision)
{
	gsl_complex sigmaKaonPropagator = 
	antiNeutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return sigmaKaonPropagator;
}


ComplexSquareMatrixGSL neutral038ScalarsPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          			   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          			   double effMassU, double effMassD, double effMassS, 
					  		          			   double k0, double k, double Gamma, double integralPrecision)
{   
	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//Set matrix with meson projection operators
	SU3NJL3DCutoffMesonProjector projector(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, integralPrecision);
	double S00 = projector.scalar00();
	double S03 = projector.scalar03();
	double S08 = projector.scalar08();
	double S33 = projector.scalar33();
	double S38 = projector.scalar38();
	double S88 = projector.scalar88();

	ComplexSquareMatrixGSL sAB(3);

	sAB.setValue(0, 0, gsl_complex_rect(2.0*S00, 0.0));
	sAB.setValue(0, 1, gsl_complex_rect(2.0*S03, 0.0));
	sAB.setValue(0, 2, gsl_complex_rect(2.0*S08, 0.0));

	sAB.setValue(1, 0, gsl_complex_rect(2.0*S03, 0.0));
	sAB.setValue(1, 1, gsl_complex_rect(2.0*S33, 0.0));
	sAB.setValue(1, 2, gsl_complex_rect(2.0*S38, 0.0));

	sAB.setValue(2, 0, gsl_complex_rect(2.0*S08, 0.0));
	sAB.setValue(2, 1, gsl_complex_rect(2.0*S38, 0.0));
	sAB.setValue(2, 2, gsl_complex_rect(2.0*S88, 0.0));

	//Set matrix with quark polarization operators
	gsl_complex scalarPolUU = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotU, effChemPotU, effMassU, effMassU, k0, k, Gamma, integralPrecision);
	gsl_complex scalarPolDD = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotD, effChemPotD, effMassD, effMassD, k0, k, Gamma, integralPrecision);
	gsl_complex scalarPolSS = 
	scalarPolarizationOperator3DCutoff(reguScheme, cutoff, Nc, T, effChemPotS, effChemPotS, effMassS, effMassS, k0, k, Gamma, integralPrecision);

	//PI00 = (2./3.)*( PIuu + PIdd + PIss)
	gsl_complex scalarPol00;
	scalarPol00 = gsl_complex_add(scalarPolUU, scalarPolDD);
    scalarPol00 = gsl_complex_add(scalarPol00, scalarPolSS);
    scalarPol00 = gsl_complex_mul_real(scalarPol00, 2.0/3.0);

    //PI03 = sqrt(2./3.)*( PIuu - PIdd )
    gsl_complex scalarPol03;
    scalarPol03 = gsl_complex_sub(scalarPolUU, scalarPolDD);
    scalarPol03 = gsl_complex_mul_real(scalarPol03, sqrt(2.0/3.0));

    //PI08 = (sqrt(2)/3.)*( PIuu + PIdd - 2 PIss )
    gsl_complex scalarPol08;
    scalarPol08 = gsl_complex_add(scalarPolUU, scalarPolDD);
    scalarPol08 = gsl_complex_add(scalarPol08, gsl_complex_mul_real(scalarPolSS, -2.0));
    scalarPol08 = gsl_complex_mul_real(scalarPol08, sqrt(2.0)/3.0);

	//PI33 = ( PIuu + PIdd )
    gsl_complex scalarPol33;
    scalarPol33 = gsl_complex_add(scalarPolUU, scalarPolDD); 

    //PI38 = (1/sqrt(3.))*( PIuu - PIdd )
    gsl_complex scalarPol38;
    scalarPol38 = gsl_complex_sub(scalarPolUU, scalarPolDD);
    scalarPol38 = gsl_complex_mul_real(scalarPol38, 1.0/sqrt(3.0));

    //PI88 = (1./3.)*( Piuu + PIdd + 4 PIss )
    gsl_complex scalarPol88;
    scalarPol88 = gsl_complex_add(scalarPolUU, scalarPolDD);
    scalarPol88 = gsl_complex_add(scalarPol88, gsl_complex_mul_real(scalarPolSS, 4.0));
    scalarPol88 = gsl_complex_mul_real(scalarPol88, 1.0/3.0);

    ComplexSquareMatrixGSL scalarPol(3);

    scalarPol.setValue(0, 0, scalarPol00);
    scalarPol.setValue(0, 1, scalarPol03);
    scalarPol.setValue(0, 2, scalarPol08);

    scalarPol.setValue(1, 0, scalarPol03);
    scalarPol.setValue(1, 1, scalarPol33);
    scalarPol.setValue(1, 2, scalarPol38);

    scalarPol.setValue(2, 0, scalarPol08);
    scalarPol.setValue(2, 1, scalarPol38);
    scalarPol.setValue(2, 2, scalarPol88);

    //Set inverse meson propagator
    ComplexSquareMatrixGSL inverseNeutral038Propagator = subtract(sAB.inverse(), scalarPol);

    //return propagator, invert the inverse propagator!
	return inverseNeutral038Propagator.inverse();
}


ComplexSquareMatrixGSL neutral038ScalarsPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          			   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          			   double effMassU, double effMassD, double effMassS, 
					  		          			   double k0, double k, double integralPrecision)
{
	ComplexSquareMatrixGSL neutral038Propagator = 
	neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision);
    
    return neutral038Propagator;
}



gsl_complex nonDiagonalMesonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		           double effChemPotU, double effChemPotD, double effChemPotS, 
					  		           double effMassU, double effMassD, double effMassS, 
					  		           double k0, double k, double Gamma, double integralPrecision,
					  		           mesonState meson)
{   
	gsl_complex mesonPropagator;

	if ( meson==pionPlus )
	{
		mesonPropagator = pionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==pionMinus )
	{
		mesonPropagator = pionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==kaonPlus )
	{
		mesonPropagator = kaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==kaonMinus )
	{
		mesonPropagator = kaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==neutralKaon )
	{
		mesonPropagator = neutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==antiNeutralKaon )
	{
		mesonPropagator = antiNeutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==sigmaPionPlus )
	{
		mesonPropagator = sigmaPionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==sigmaPionMinus )
	{
		mesonPropagator = sigmaPionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==sigmaKaonPlus )
	{
		mesonPropagator = sigmaKaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==sigmaKaonMinus )
	{
		mesonPropagator = sigmaKaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==neutralSigmaKaon )
	{
		mesonPropagator = neutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if ( meson==antiNeutralSigmaKaon )
	{
		mesonPropagator = antiNeutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}

    return mesonPropagator;
}


gsl_complex nonDiagonalMesonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		           double effChemPotU, double effChemPotD, double effChemPotS, 
					  		           double effMassU, double effMassD, double effMassS, 
					  		           double k0, double k, double integralPrecision,
					  		           mesonState meson)
{
	gsl_complex mesonPropagator = 
	nonDiagonalMesonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision, meson);
    
    return mesonPropagator;
}


ComplexSquareMatrixGSL diagonalMesonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          		   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          		   double effMassU, double effMassD, double effMassS, 
					  		          		   double k0, double k, double Gamma, double integralPrecision,
					  		          		   mesonState meson)
{
	ComplexSquareMatrixGSL mesonPropagator(3);

	if( meson==diagonalPseudoscalars )
	{
		mesonPropagator = 
		neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}
	else if( meson==diagonalScalars )
	{
		mesonPropagator = 
		neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, Gamma, integralPrecision);
	}

	return mesonPropagator;
}


ComplexSquareMatrixGSL diagonalMesonPropagator(SU3NJL3DCutoffParameters parametersNJL, double T, 
					  		          		   double effChemPotU, double effChemPotD, double effChemPotS, 
					  		          		   double effMassU, double effMassD, double effMassS, 
					  		          		   double k0, double k, double integralPrecision,
					  		          		   mesonState meson)
{
	ComplexSquareMatrixGSL mesonPropagator = 
	diagonalMesonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, k0, k, 0.0, integralPrecision, meson);

	return mesonPropagator;
}


int SU3NJL3DCutoffMesonMassEquations(const gsl_vector *x, void *auxiliar, gsl_vector *f)
{
	//define variables
    double mesonMass = gsl_vector_get(x,0);
    double mesonWidth = gsl_vector_get(x,1);

    //define parameters
    SU3NJL3DCutoffMeson meson(auxiliar);

    //calculate inverse meson propagators and set the real and imaginary parts to zero
   	double k0 = mesonMass;
	double k = 0.0;
	double gamma = mesonWidth;
	double f0, f1;
	if ( meson.getMesonState()!=diagonalPseudoscalars && meson.getMesonState()!=diagonalScalars )
	{
	    gsl_complex inverseMesonPropagator = meson.calculateInverseNonDiagonalPropagator(k0, k, gamma);

	    f0 = GSL_REAL( inverseMesonPropagator );
	    f1 = GSL_IMAG( inverseMesonPropagator );
	}
	else
	{
	    vector<gsl_complex> eigenvalue = meson.calculateInverseDiagonalPropagatorEigenvalues(k0, k, gamma);

	    f0 = GSL_REAL( eigenvalue[0] );
	    f1 = GSL_IMAG( eigenvalue[0] );
	}
	gsl_vector_set (f, 0, f0);
	gsl_vector_set (f, 1, f1);

	return GSL_SUCCESS;
}






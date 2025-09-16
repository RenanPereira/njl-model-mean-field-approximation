#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.h"
#include "njl_model/n_fermion_line_integrals/one_fermion_line_integral_3d_cutoff.h"
#include "gsl_wrapper/root_solver_gsl.h"


//Zero variables necessary in this file: ICS=IntegratedCrossSections
const double ICS_ZERO = 1E-15;
const double ICS_ZERO_MASS_DIFFERENCE = 1E-8;


//In this file we have functions to calculate quark-quark and quark-antiquark integrated cross sections
//quark-quark: 
//UD->UD, DU->DU, US->US, SU->SU, DS->DS, SD->SD, UU->UU, DD->DD, SS->SS
//quark-antiquark: 
//UDBar->UDBar, USBar->USBar, DSBar->DSBar,
//DUBar->DUBar, SUBar->SUBar, SDBar->SDBar,
//UUBar->UUBar, UUBar->DDBar, UUBar->SSBar,
//DDBar->UUBar, DDBar->DDBar, DDBar->SSBar,
//SSBar->UUBar, SSBar->DDBar, SSBar->SSBar
//UBarUBar->UBarUBar, UBarDBar->UBarDBar, UBarSBar->UBarSBar, SBarSBar->SBarSBar


////////////////////////////////////////////////////////////////////////////////////////
//Fermi-Dirac function integration: "number density" (not to be confused with quark density)
double fermiDiracIntegralIntegrand(double p, void *parameters)
{   
	OneFermionLine3DCutoffIntegrand aux(parameters);
	double T = aux.getTemperature();
	double effCP = aux.getEffectiveChemicalPotential();
	double M = aux.getEffectiveMass();

	double E = sqrt( pow(p,2) + pow(M,2) );

	double fermiDiracIntegrand = pow(p,2)*fermiDistribution(T, E-effCP);

    return fermiDiracIntegrand;
}


double fermiDiracIntegral(double Nc, double T, double effCPot, double effMass, double cutoff, double integralPrecision)
{	
	OneFermionLine3DCutoffIntegrand aux("fermiDiracIntegral", T, effCPot, effMass);

    int integrationWorkspace = 1000;
    Integration1DimGSLQAGS fermiDiracIntegral_dp(0.0, cutoff, &aux, fermiDiracIntegralIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
    double fermiDirac = ( Nc/(M_PI*M_PI) )*fermiDiracIntegral_dp.evaluate();

    return fermiDirac;
}


double fermiDiracIntegral(double Nc, double T, double effCPot, double effMass, double integralPrecision)
{	
	OneFermionLine3DCutoffIntegrand aux("fermiDiracIntegral", T, effCPot, effMass);

    int integrationWorkspace = 1000;
    Integration1DimGSLQAGIU fermiDiracIntegral_dp(0.0, &aux, fermiDiracIntegralIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
    double fermiDirac = ( Nc/(M_PI*M_PI) )*fermiDiracIntegral_dp.evaluate();

    return fermiDirac;
}


////////////////////////////////////////////////////////////////////////////////////////
//12->34 average transition rate W: integrating up to the cutoff with change of variables


double dsdEdepsilon_epsilonMinPlus(double M1, double M2, double s, double E)
{	
	double epsilonMinPlus = 0.0;

	if ( fabs(M1-M2)>ICS_ZERO_MASS_DIFFERENCE )
	{	
		double sqrtArg = ( pow(M1, 4) + pow(pow(M2,2)-s, 2) - 2.0*pow(M1, 2)*( pow(M2, 2) + s ) )*( pow(E,2) - s );

		epsilonMinPlus = ( -pow(M1,2)*E + pow(M2,2)*E + sqrt( fabs(sqrtArg) ) )/s;

		if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
		{ 
			cout << "Argument of the square root in dsdEdepsilon_epsilonMinPlus is negative and larger then variable ICS_ZERO!\n"; 
		}
	}
	else
	{
		double sqrtArg = ( s - 4.0*pow(M1,2) )*( pow(E, 2) - s );

		epsilonMinPlus = sqrt( fabs(sqrtArg) )/sqrt(s);

		if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
		{ 
			cout << "Argument of the square root in dsdEdepsilon_epsilonMinPlus is negative and larger then variable ICS_ZERO!\n"; 
		}
	}

	return epsilonMinPlus;
}


double dsdEdepsilon_epsilonMinMinus(double M1, double M2, double s, double E)
{	
	double epsilonMinMinus = 0.0;

	if ( fabs(M1-M2)>ICS_ZERO_MASS_DIFFERENCE )
	{	
		double sqrtArg = ( pow(M1, 4) + pow(pow(M2,2)-s, 2) - 2.0*pow(M1, 2)*( pow(M2, 2) + s ) )*( pow(E,2) - s );

		epsilonMinMinus = ( -pow(M1,2)*E + pow(M2,2)*E - sqrt( fabs(sqrtArg) ) )/s;

		if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
		{ 
			cout << "Argument of the square root in dsdEdepsilonepsilonminminus is negative and larger then variable ICS_ZERO!\n"; 
		}
	}
	else
	{
		double sqrtArg = ( s - 4.0*pow(M1,2) )*( pow(E, 2) - s );

		epsilonMinMinus = -sqrt( fabs(sqrtArg) )/sqrt(s);

		if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
		{ 
			cout << "Argument of the square root in dsdEdepsilonepsilonminminus is negative and larger then variable ICS_ZERO!\n"; 
		}
	}

	return epsilonMinMinus;
}


double dsdEdepsilon_epsilonLambdaM1Plus(double cutoff, double M1, double E)
{
	double epsilonLambdaM1Plus = E + 2.0*sqrt( pow(cutoff, 2) + pow(M1, 2) );

	return epsilonLambdaM1Plus;
}


double dsdEdepsilon_epsilonLambdaM1Minus(double cutoff, double M1, double E)
{
	double epsilonLambdaM1Minus = E - 2.0*sqrt( pow(cutoff, 2) + pow(M1, 2) );

	return epsilonLambdaM1Minus;
}


double dsdEdepsilon_epsilonLambdaM2Plus(double cutoff, double M2, double E)
{
	double epsilonLambdaM2Plus = -E + 2.0*sqrt( pow(cutoff, 2) + pow(M2, 2) );

	return epsilonLambdaM2Plus;
}


double dsdEdepsilon_epsilonLambdaM2Minus(double cutoff, double M2, double E)
{
	double epsilonLambdaM2Minus = -E - 2.0*sqrt( pow(cutoff, 2) + pow(M2, 2) );

	return epsilonLambdaM2Minus;
}


double dsdEdepsilon_gamma1E(double cutoff, double M1, double M2, double s)
{	
	double sqrtArg = ( pow(M1-M2,2) - s )*( pow(M1+M2,2) - s )*pow(cutoff,2);

	double gamma1E = ( -sqrt( fabs(sqrtArg) ) + ( pow(M2,2) - pow(M1,2) + s )*sqrt( pow(cutoff,2) + pow(M2,2) ) )/( 2.0*pow(M2,2) );

	if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
	{ 
		cout << "Argument of the square root in dsdEdepsilon_gamma1E is negative and larger then variable ICS_ZERO!\n"; 
	}

	return gamma1E;
}


double dsdEdepsilon_beta1E(double cutoff, double M1, double M2, double s)
{	
	double sqrtArg = ( pow(M1-M2,2) - s )*( pow(M1+M2,2) - s )*pow(cutoff,2);

	double beta1E = ( sqrt( fabs(sqrtArg) ) + ( pow(M2,2) - pow(M1,2) + s )*sqrt( pow(cutoff,2) + pow(M2,2) ) )/( 2.0*pow(M2,2) );

	if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
	{ 
		cout << "Argument of the square root in dsdEdepsilon_beta1E is negative and larger then variable ICS_ZERO!\n"; 
	}

	return beta1E;
}


double dsdEdepsilon_gamma1epsilon(double cutoff, double M1, double M2, double s)
{	
	double sqrtArg = ( pow(M1-M2,2) - s )*( pow(M1+M2,2) - s )*pow(cutoff,2);

	double gamma1epsilon = ( sqrt( fabs(sqrtArg) ) + ( pow(M1,2) + 3.0*pow(M2,2) - s )*sqrt( pow(cutoff,2) + pow(M2,2) ) )/( 2.0*pow(M2,2) );

	if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
	{ 
		cout << "Argument of the square root in dsdEdepsilon_gamma1epsilon is negative and larger then variable ICS_ZERO!\n"; 
	}

	return gamma1epsilon;
}


double dsdEdepsilon_beta1epsilon(double cutoff, double M1, double M2, double s)
{	
	double sqrtArg = ( pow(M1-M2,2) - s )*( pow(M1+M2,2) - s )*pow(cutoff,2);

	double beta1epsilon = ( -sqrt( fabs(sqrtArg) ) + ( pow(M1,2) + 3.0*pow(M2,2) - s )*sqrt( pow(cutoff,2) + pow(M2,2) ) )/( 2.0*pow(M2,2) );

	if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
	{ 
		cout << "Argument of the square root in dsdEdepsilonbeta1epsilon is negative and larger then variable ICS_ZERO!\n"; 
	}

	return beta1epsilon;
}


double dsdEdepsilon_gamma2E(double cutoff, double M1, double M2, double s)
{	
	double sqrtArg = ( pow(M1-M2,2) - s )*( pow(M1+M2,2) - s )*pow(cutoff,2);

	double gamma2E = ( -sqrt( fabs(sqrtArg) ) + ( pow(M1,2) - pow(M2,2) + s )*sqrt( pow(cutoff,2) + pow(M1,2) ) )/( 2.0*pow(M1,2) );

	if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
	{ 
		cout << "Argument of the square root in dsdEdepsilon_gamma2E is negative and larger then variable ICS_ZERO!\n"; 
	}

	return gamma2E;
}


double dsdEdepsilon_beta2E(double cutoff, double M1, double M2, double s)
{	
	double sqrtArg = ( pow(M1-M2,2) - s )*( pow(M1+M2,2) - s )*pow(cutoff,2);

	double beta2E = ( sqrt( fabs(sqrtArg) ) + ( pow(M1,2) - pow(M2,2) + s )*sqrt( pow(cutoff,2) + pow(M1,2) ) )/( 2.0*pow(M1,2) );

	if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
	{ 
		cout << "Argument of the square root in dsdEdepsilon_beta2E is negative and larger then variable ICS_ZERO!\n"; 
	}

	return beta2E;
}


double dsdEdepsilon_gamma2epsilon(double cutoff, double M1, double M2, double s)
{	
	double sqrtArg = ( pow(M1-M2,2) - s )*( pow(M1+M2,2) - s )*pow(cutoff,2);

	double gamma2epsilon = ( -sqrt( fabs(sqrtArg) ) - ( 3.0*pow(M1,2) + pow(M2,2) - s )*sqrt( pow(cutoff,2) + pow(M1,2) ) )/( 2.0*pow(M1,2) );

	if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
	{ 
		cout << "Argument of the square root in dsdEdepsilon_gamma2epsilon is negative and larger then variable ICS_ZERO!\n"; 
	}

	return gamma2epsilon;
}


double dsdEdepsilon_beta2epsilon(double cutoff, double M1, double M2, double s)
{	
	double sqrtArg = ( pow(M1-M2,2) - s )*( pow(M1+M2,2) - s )*pow(cutoff,2);

	double beta2epsilon = ( sqrt( fabs(sqrtArg) ) - ( 3.0*pow(M1,2) + pow(M2,2) - s )*sqrt( pow(cutoff,2) + pow(M1,2) ) )/( 2.0*pow(M1,2) );

	if ( sqrtArg<0 && fabs(sqrtArg)>ICS_ZERO )
	{ 
		cout << "Argument of the square root in dsdEdepsilon_beta2epsilon is negative and larger then variable ICS_ZERO!\n"; 
	}

	return beta2epsilon;
}


double dsdEdepsilon_alpha1E(double s)
{
	double alpha1E = sqrt( s );

	return alpha1E;
}


double dsdEdepsilon_alpha1epsilon(double M1, double M2, double s)
{
	double alpha1epsilon = ( pow(M2,2) - pow(M1,2) )/sqrt( s );

	return alpha1epsilon;
}


double dsdEdepsilon_ELambdaSwitchE(double cutoff, double M1, double M2)
{
	double ELambdaSwitchE = sqrt( pow(cutoff,2) + pow(M1,2) ) + sqrt( pow(cutoff,2) + pow(M2,2) );

	return ELambdaSwitchE;
}


double dsdEdepsilon_ELambdaswitchepsilon(double cutoff, double M1, double M2)
{
	double ELambdaswitchepsilon = -sqrt( pow(cutoff,2) + pow(M1,2) ) + sqrt( pow(cutoff,2) + pow(M2,2) );

	return ELambdaswitchepsilon;
}


double dsdEdepsilon_epsilonMax(double cutoff, double M1, double M2, double s, double E)
{	
	double epsilonMax = 0.0;

	if ( E > dsdEdepsilon_gamma1E(cutoff, M1, M2, s) )
	{
		epsilonMax = dsdEdepsilon_epsilonLambdaM2Plus(cutoff, M2, E);
	}
	else
	{
		epsilonMax = dsdEdepsilon_epsilonMinPlus(M1, M2, s, E);
	}

	return epsilonMax;
}


double dsdEdepsilon_epsilonMin(double cutoff, double M1, double M2, double s, double E)
{	
	double epsilonMin = 0.0;

	if ( E > dsdEdepsilon_gamma2E(cutoff, M1, M2, s) )
	{
		epsilonMin = dsdEdepsilon_epsilonLambdaM1Minus(cutoff, M1, E);
	}
	else
	{
		epsilonMin = dsdEdepsilon_epsilonMinMinus(M1, M2, s, E);
	}

	return epsilonMin;
}


double dsdEdepsilon_EMinM1LargerM2(double cutoff, double M1, double M2, double s)
{
	double EMinM1LargerM2 = 0.0;

	double alpha1E = dsdEdepsilon_alpha1E(s);

	if ( dsdEdepsilon_epsilonLambdaM1Minus(cutoff, M1, alpha1E)>dsdEdepsilon_alpha1epsilon(M1, M2, s) )
	{
		EMinM1LargerM2 = dsdEdepsilon_gamma2E(cutoff, M1, M2, s);
	}
	else
	{
		EMinM1LargerM2 = alpha1E;
	}

	return EMinM1LargerM2;
}


double dsdEdepsilon_EMinM2LargerM1(double cutoff, double M1, double M2, double s)
{
	double EMinM2LargerM1 = 0.0;

	double alpha1E = dsdEdepsilon_alpha1E(s);

	if ( dsdEdepsilon_epsilonLambdaM2Plus(cutoff, M2, alpha1E)<dsdEdepsilon_alpha1epsilon(M1, M2, s) )
	{
		EMinM2LargerM1 = dsdEdepsilon_gamma1E(cutoff, M1, M2, s);
	}
	else
	{
		EMinM2LargerM1 = alpha1E;
	}

	return EMinM2LargerM1;
}


double dsdEdepsilon_EMin(double cutoff, double M1, double M2, double s)
{
	double EMin = 0.0;

	if ( s>pow(dsdEdepsilon_ELambdaSwitchE(cutoff, M1, M2),2) )
	{
		EMin = 0.0;
	}
	else
	{
		if ( fabs(M1-M2)>ICS_ZERO_MASS_DIFFERENCE )
		{
			if ( M1>M2 )
			{
				EMin = dsdEdepsilon_EMinM1LargerM2(cutoff, M1, M2, s);
			}
			else
			{
				EMin = dsdEdepsilon_EMinM2LargerM1(cutoff, M1, M2, s);
			}
		}
		else
		{
			EMin = dsdEdepsilon_alpha1E(s);
		}
	}

	return EMin;
}


double dsdEdepsilon_EMaxM1LargerM2(double cutoff, double M1, double M2, double s)
{
	double EMaxM1LargerM2 = 0.0;

	double ELambdaSwitchE = dsdEdepsilon_ELambdaSwitchE(cutoff, M1, M2);
	double beta2E = dsdEdepsilon_beta2E(cutoff, M1, M2, s);

	if ( ELambdaSwitchE<beta2E )
	{
		EMaxM1LargerM2 = ELambdaSwitchE;
	}
	else
	{
		EMaxM1LargerM2 = beta2E;
	}

	return EMaxM1LargerM2;
}


double dsdEdepsilon_EMaxM2LargerM1(double cutoff, double M1, double M2, double s)
{
	double EMaxM2LargerM1 = 0.0;

	double ELambdaSwitchE = dsdEdepsilon_ELambdaSwitchE(cutoff, M1, M2);
	double beta1E = dsdEdepsilon_beta1E(cutoff, M1, M2, s);

	if ( ELambdaSwitchE<beta1E )
	{
		EMaxM2LargerM1 = ELambdaSwitchE;
	}
	else
	{
		EMaxM2LargerM1 = beta1E;
	}

	return EMaxM2LargerM1;
}


double dsdEdepsilon_EMax(double cutoff, double M1, double M2, double s)
{
	double EMax = 0.0;

	if ( s>pow(dsdEdepsilon_ELambdaSwitchE(cutoff, M1, M2),2) )
	{
		EMax = 0.0;
	}
	else
	{
		if ( fabs(M1-M2)>ICS_ZERO_MASS_DIFFERENCE )
		{
			if ( M1>M2 )
			{
				EMax = dsdEdepsilon_EMaxM1LargerM2(cutoff, M1, M2, s);
			}
			else
			{
				EMax = dsdEdepsilon_EMaxM2LargerM1(cutoff, M1, M2, s);
			}
		}
		else
		{
			EMax = dsdEdepsilon_ELambdaSwitchE(cutoff, M1, M2);
		}
	}

	return EMax;
}


double integratedCrossSectionCOVVolumeIntegrand_dsdE(double E, void *parameters)
{
	IntegratedCrossSectionIntegrand aux(parameters);
    double cutoff = aux.getParametersNJL().getThreeMomentumCutoff();
	double M1 = aux.getUpQuarkEffectiveMass();
	double M2 = aux.getDownQuarkEffectiveMass();
	double s = aux.getCenterOfMassEnergy();

    double volumeIntegrand = dsdEdepsilon_epsilonMax(cutoff, M1, M2, s, E) - dsdEdepsilon_epsilonMin(cutoff, M1, M2, s, E);

    return volumeIntegrand;
}


double integratedCrossSectionCOVVolumeIntegrand_ds(double s, void *parameters)
{
	IntegratedCrossSectionIntegrand aux(parameters);
    double cutoff = aux.getParametersNJL().getThreeMomentumCutoff();
	double M1 = aux.getUpQuarkEffectiveMass();
	double M2 = aux.getDownQuarkEffectiveMass();
	double integralPrecision_dsdE = aux.getIntegratedCrossSectionIntegralPrecision_dXdY();

	aux.setIntegralID("integratedCrossSectionCOVVolume_dsdE");
	aux.setCenterOfMassEnergy(s);

    double EMin = dsdEdepsilon_EMin(cutoff, M1, M2, s);
    double EMax = dsdEdepsilon_EMax(cutoff, M1, M2, s);
    
    int integrationWorkspace = 1000;
    Integration1DimGSLQAG integratedCrossSectionCOVVolume_dsdE(EMin, EMax, &aux, integratedCrossSectionCOVVolumeIntegrand_dsdE, integralPrecision_dsdE, integralPrecision_dsdE, integrationWorkspace, 1);
    double volumeIntegrand = integratedCrossSectionCOVVolume_dsdE.evaluate();

	return volumeIntegrand;
}


double integratedCrossSectionCOVVolume(double cutoff, double M1, double M2, double integralPrecision_dsdE, double integralPrecision_ds)
{	
	SU3NJL3DCutoffParameters parametersNJL(cutoff);
	IntegratedCrossSectionIntegrand aux(
		"integratedCrossSectionCOVVolumeIntegral", parametersNJL, 0.0, 
        0.0, 0.0, 0.0, 
        M1, M2, 0.0, 
        0.0, UDUD, 
        false, 0.0,
        integralPrecision_dsdE
	);

	double sMin = pow(M1+M2,2);
    double sMax = pow(sqrt( pow(cutoff,2) + pow(M1,2) ) + sqrt( pow(cutoff,2) + pow(M2,2) ), 2);

    int integrationWorkspace = 1000;
    Integration1DimGSLQAG integratedCrossSectionCOVVolumeIntegral(sMin, sMax, &aux, integratedCrossSectionCOVVolumeIntegrand_ds, integralPrecision_ds, integralPrecision_ds, integrationWorkspace, 1);
    double volumeIntegrand = (1.0/4.0)*integratedCrossSectionCOVVolumeIntegral.evaluate();

    return volumeIntegrand;
}


double integratedCrossSectionOGVolumeIntegrand_dp1dp2dtheta(double theta, void *parameters)
{
	IntegratedCrossSectionIntegrand aux(parameters);
	double M1 = aux.getUpQuarkEffectiveMass();
	double M2 = aux.getDownQuarkEffectiveMass();
	double p1 = aux.getMomentumParticle1();
	double p2 = aux.getMomentumParticle2();

	double E1 = sqrt( pow(M1,2) + pow(p1,2) );
	double E2 = sqrt( pow(M2,2) + pow(p2,2) );

    double volumeIntegrand = ( pow(p1,2)*pow(p2,2)*sin(theta) )/( E1*E2 );

    return volumeIntegrand;
}


double integratedCrossSectionOGVolumeIntegrand_dp1dp2(double p2, void *parameters)
{
	IntegratedCrossSectionIntegrand aux(parameters);
	double integralPrecision = aux.getIntegratedCrossSectionIntegralPrecision_dXdYdZ();

	aux.setIntegralID("integratedCrossSectionIntegral_dp1dp2dtheta");
	aux.setMomentumParticle2(p2);
	
	int integrationWorkspace = 1000;
    Integration1DimGSLQAGS integratedCrossSectionIntegral_dp1dp2dtheta(0, M_PI, &aux, integratedCrossSectionOGVolumeIntegrand_dp1dp2dtheta, integralPrecision, integralPrecision, integrationWorkspace);
    double volumeIntegrand = integratedCrossSectionIntegral_dp1dp2dtheta.evaluate();
	
	return volumeIntegrand;

	/*
	IntegratedCrossSectionIntegrand aux(parameters);

	double M1 = aux.getUpQuarkEffectiveMass();
	double M2 = aux.getDownQuarkEffectiveMass();
	double p1 = aux.getMomentumParticle1();

	double E1 = sqrt( pow(M1,2) + pow(p1,2) );
	double E2 = sqrt( pow(M2,2) + pow(p2,2) );

    double volumeIntegrand = 2.0*( pow(p1,2)*pow(p2,2) )/( E1*E2 );

    return volumeIntegrand;
    */
}


double integratedCrossSectionOGVolumeIntegrand_dp1(double p1, void *parameters)
{
	IntegratedCrossSectionIntegrand aux(parameters);
	double cutoff = aux.getParametersNJL().getThreeMomentumCutoff();
	double integralPrecision = aux.getIntegratedCrossSectionIntegralPrecision_dXdY();

	aux.setIntegralID("integratedCrossSectionIntegral_dp1dp2");
	aux.setMomentumParticle1(p1);

	int integrationWorkspace = 1000;
    Integration1DimGSLQAGS integratedCrossSectionIntegral_dp1dp2(0, cutoff, &aux, integratedCrossSectionOGVolumeIntegrand_dp1dp2, integralPrecision, integralPrecision, integrationWorkspace);
    double volumeIntegrand = integratedCrossSectionIntegral_dp1dp2.evaluate();

    return volumeIntegrand;
}


double integratedCrossSectionOGVolume(double cutoff, double M1, double M2, double integralPrecision_dp1dp2dtheta, double integralPrecision_dp1dp2, double integralPrecision_dp1)
{
	SU3NJL3DCutoffParameters parametersNJL(cutoff);
	IntegratedCrossSectionIntegrand aux(
		"integratedCrossSectionOGVolumeIntegral", parametersNJL, 0.0, 
        0.0, 0.0, 0.0, 
        M1, M2, 0.0, 
        0.0, UDUD, 
        false, 0.0,
        integralPrecision_dp1dp2dtheta, integralPrecision_dp1dp2
	);

    int integrationWorkspace = 1000;
    Integration1DimGSLQAGS integratedCrossSectionOGVolumeIntegral(0, cutoff, &aux, integratedCrossSectionOGVolumeIntegrand_dp1, integralPrecision_dp1, integralPrecision_dp1, integrationWorkspace);
    double volumeIntegrand = integratedCrossSectionOGVolumeIntegral.evaluate();

    return volumeIntegrand;
}


double nEta(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double s, double E)
{   
    double C = dsdEdepsilon_epsilonMin(cutoff, M1, M2, s, E);
    double D = dsdEdepsilon_epsilonMax(cutoff, M1, M2, s, E);
    double eta = -1.0;

/*
    double gpluseta_aux;
    gpluseta_aux = 2*bosedist(T, E+eta*(effCP1+effCP2))*( 
													+ T*log( 1 + exp( 0.5*(C - E - 2*eta*effCP1)/T ) )
													- T*log( 1 + exp( 0.5*(D - E - 2*eta*effCP1)/T ) )
													+ T*log( 1 + exp( 0.5*(D + E + 2*eta*effCP2)/T ) )
													- T*log( 1 + exp( 0.5*(C + E + 2*eta*effCP2)/T ) )
    											);
*/

    double argC1 = 0.5*(C - E - 2.0*eta*effCP1);
    double argD1 = 0.5*(D - E - 2.0*eta*effCP1);
    double argD2 = 0.5*(D + E + 2.0*eta*effCP2);
    double argC2 = 0.5*(C + E + 2.0*eta*effCP2);

    double nEtaAux = 0.0;
    {
        if ( isinf( exp(argC1/T) )==0 )
        { 
            nEtaAux = nEtaAux + T*log( 1.0 + exp(argC1/T) ); 
        }
        else if( argC1>0 )
        { 
            nEtaAux = nEtaAux + puiseuxExpansionTln1plusExpArgOverT(T, argC1); 
        }

        if ( isinf( exp(argD1/T) )==0 )
        { 
            nEtaAux = nEtaAux - T*log( 1.0 + exp(argD1/T) ); 
        }
        else if( argD1>0 )
        { 
            nEtaAux = nEtaAux - puiseuxExpansionTln1plusExpArgOverT(T, argD1); 
        }

        if ( isinf( exp(argD2/T) )==0 )
        { 
            nEtaAux = nEtaAux + T*log( 1.0 + exp(argD2/T) ); 
        }
        else if( argD2>0 )
        { 
            nEtaAux = nEtaAux + puiseuxExpansionTln1plusExpArgOverT(T, argD2); 
        }

        if ( isinf( exp(argC2/T) )==0 )
        { 
            nEtaAux = nEtaAux - T*log( 1.0 + exp(argC2/T) ); 
        }
        else if( argC2>0 )
        { 
            nEtaAux = nEtaAux - puiseuxExpansionTln1plusExpArgOverT(T, argC2); 
        }
    }

    nEtaAux = 2.0*boseDistribution(T, E+eta*(effCP1+effCP2))*nEtaAux;

    return nEtaAux;
}


double integratedCrossSectionIntegrand_dsdE(double E, void *parameters)
{	
	IntegratedCrossSectionIntegrand aux(parameters);
	double cutoff = aux.getParametersNJL().getThreeMomentumCutoff();
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	double effCPD = aux.getDownQuarkEffectiveChemicalPotential();
	double effCPS = aux.getStrangeQuarkEffectiveChemicalPotential();
	double effMU = aux.getUpQuarkEffectiveMass();
	double effMD = aux.getDownQuarkEffectiveMass();
	double effMS = aux.getStrangeQuarkEffectiveMass();
	double s = aux.getCenterOfMassEnergy();
	scatteringProcess process = aux.getProcess();

	//set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMU, effMD, effMS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effCPU, effCPD, effCPS, process, cP1, cP2);

    double integratedCrossSection_dsdE = nEta(T, cP1, cP2, cutoff, m1, m2, s, E);

    return integratedCrossSection_dsdE;
}


double integratedCrossSectionIntegrand_ds(double s, void *parameters)
{
	IntegratedCrossSectionIntegrand aux(parameters);
	double cutoff = aux.getParametersNJL().getThreeMomentumCutoff();
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	double effCPD = aux.getDownQuarkEffectiveChemicalPotential();
	double effCPS = aux.getStrangeQuarkEffectiveChemicalPotential();
	double effMU = aux.getUpQuarkEffectiveMass();
	double effMD = aux.getDownQuarkEffectiveMass();
	double effMS = aux.getStrangeQuarkEffectiveMass();
	scatteringProcess process = aux.getProcess();
	double integralPrecision_dsdE = aux.getIntegratedCrossSectionIntegralPrecision_dXdY();

	//set center of mass energy
	aux.setCenterOfMassEnergy(s);

    //set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMU, effMD, effMS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effCPU, effCPD, effCPS, process, cP1, cP2);
	
	//calculate cross section	
	double crossSection = 1.0;
	crossSection = crossSectionProcess12To34(
		aux.getParametersNJL(), T, 
        effCPU, effCPD, effCPS, 
        effMU, effMD, effMS, 
        s, aux.getPropagatorIntegralPrecision(), process, 
        aux.getLargeAngleScatteringContribution(), aux.getCrossSectionIntegralPrecision()
	);

    double EMin = dsdEdepsilon_EMin(cutoff, m1, m2, s);
    double EMax = dsdEdepsilon_EMax(cutoff, m1, m2, s);

    double g1 = dsdEdepsilon_gamma1E(cutoff, m1, m2, s);
    double g2 = dsdEdepsilon_gamma2E(cutoff, m1, m2, s);


	int integrationWorkspace = 1000;
    double integratedCrossSection_ds = 0.0;

	if ( m2>m1 && fabs(m1-m2)>ICS_ZERO_MASS_DIFFERENCE ) //in this case g2>=g1
    {   
    	aux.setIntegralID("integratedCrossSectionIntegral_dsdE_EMin_g1");
    	Integration1DimGSLQAGS integratedCrossSectionIntegral_dsdE_EMin_g1(EMin, g1, &aux, integratedCrossSectionIntegrand_dsdE, integralPrecision_dsdE, integralPrecision_dsdE, integrationWorkspace);
    	if ( fabs(g1-EMin)>ICS_ZERO )
    	{
    		integratedCrossSection_ds = integratedCrossSection_ds + integratedCrossSectionIntegral_dsdE_EMin_g1.evaluate();
    	}

    	aux.setIntegralID("integratedCrossSectionIntegral_dsdE_g1_g2");
    	Integration1DimGSLQAGS integratedCrossSectionIntegral_dsdE_g1_g2(g1, g2, &aux, integratedCrossSectionIntegrand_dsdE, integralPrecision_dsdE, integralPrecision_dsdE, integrationWorkspace);
    	if ( fabs(g2-g1)>ICS_ZERO )
    	{
    		integratedCrossSection_ds = integratedCrossSection_ds + integratedCrossSectionIntegral_dsdE_g1_g2.evaluate();
    	}

    	aux.setIntegralID("integratedCrossSectionIntegral_dsdE_g2_EMax");
    	Integration1DimGSLQAGS integratedCrossSectionIntegral_dsdE_g2_EMax(g2, EMax, &aux, integratedCrossSectionIntegrand_dsdE, integralPrecision_dsdE, integralPrecision_dsdE, integrationWorkspace);
    	if ( fabs(EMax-g2)>ICS_ZERO )
    	{
    		integratedCrossSection_ds = integratedCrossSection_ds + integratedCrossSectionIntegral_dsdE_g2_EMax.evaluate();
    	}
    }
    else if( m1>m2 && fabs(m1-m2)>ICS_ZERO_MASS_DIFFERENCE )  //in this case g1>=g2
    {
    	aux.setIntegralID("integratedCrossSectionIntegral_dsdE_EMin_g2");
    	Integration1DimGSLQAGS integratedCrossSectionIntegral_dsdE_EMin_g2(EMin, g2, &aux, integratedCrossSectionIntegrand_dsdE, integralPrecision_dsdE, integralPrecision_dsdE, integrationWorkspace);
    	if ( fabs(g2-EMin)>ICS_ZERO )
    	{
    		integratedCrossSection_ds = integratedCrossSection_ds + integratedCrossSectionIntegral_dsdE_EMin_g2.evaluate();
    	}

    	aux.setIntegralID("integratedCrossSectionIntegral_dsdE_g2_g1");
    	Integration1DimGSLQAGS integratedCrossSectionIntegral_dsdE_g2_g1(g2, g1, &aux, integratedCrossSectionIntegrand_dsdE, integralPrecision_dsdE, integralPrecision_dsdE, integrationWorkspace);
    	if ( fabs(g1-g2)>ICS_ZERO )
    	{
    		integratedCrossSection_ds = integratedCrossSection_ds + integratedCrossSectionIntegral_dsdE_g2_g1.evaluate();
    	}

    	aux.setIntegralID("integratedCrossSectionIntegral_dsdE_g1_EMax");
    	Integration1DimGSLQAGS integratedCrossSectionIntegral_dsdE_g1_EMax(g1, EMax, &aux, integratedCrossSectionIntegrand_dsdE, integralPrecision_dsdE, integralPrecision_dsdE, integrationWorkspace);
    	if ( fabs(EMax-g1)>ICS_ZERO )
    	{
    		integratedCrossSection_ds = integratedCrossSection_ds + integratedCrossSectionIntegral_dsdE_g1_EMax.evaluate();
    	}
    }
    else //M1==M2, in this case g2==g1
    {
    	aux.setIntegralID("integratedCrossSectionIntegral_dsdE_EMin_g1");
    	Integration1DimGSLQAGS integratedCrossSectionIntegral_dsdE_EMin_g1(EMin, g1, &aux, integratedCrossSectionIntegrand_dsdE, integralPrecision_dsdE, integralPrecision_dsdE, integrationWorkspace);
    	if ( fabs(g1-EMin)>ICS_ZERO )
    	{
    		integratedCrossSection_ds = integratedCrossSection_ds + integratedCrossSectionIntegral_dsdE_EMin_g1.evaluate();
    	}
    	
    	aux.setIntegralID("integratedCrossSectionIntegral_dsdE_g1_EMax");
    	Integration1DimGSLQAGS integratedCrossSectionIntegral_dsdE_g1_EMax(g1, EMax, &aux, integratedCrossSectionIntegrand_dsdE, integralPrecision_dsdE, integralPrecision_dsdE, integrationWorkspace);
    	if ( fabs(EMax-g1)>ICS_ZERO )
    	{
    		integratedCrossSection_ds = integratedCrossSection_ds + integratedCrossSectionIntegral_dsdE_g1_EMax.evaluate();
    	}
    }

    integratedCrossSection_ds = ( sqrt(s)*momentumCM(s, m1, m2)*crossSection )*integratedCrossSection_ds; 

    cout << toString(process) << "\t" 
         << "T = " << T << "\t" 
         << "cPU = " << effCPU << "\t" 
    	 << "s = " << s << "\t" <<  "integratedCrossSection_ds = " << integratedCrossSection_ds << "\n";

    return integratedCrossSection_ds;
}


double integratedCrossSectionCOVNormalizedIntegrand_dx(double x, void *parameters)
{	
	IntegratedCrossSectionIntegrand aux(parameters);
	double normalization = aux.getNormalizationRiemannSum_ds();

	//change of variables
	double s = ( 1.0 - x )/x;

    double integratedCrossSection_dx = integratedCrossSectionIntegrand_ds(s, parameters);

    integratedCrossSection_dx = integratedCrossSection_dx/pow(x,2);

    integratedCrossSection_dx = integratedCrossSection_dx/normalization;

    return integratedCrossSection_dx;
}


double integratedCrossSectionProcess12To34(
	SU3NJL3DCutoffParameters parametersNJL, double T, 
	double effChemPotU, double effChemPotD, double effChemPotS, 
	double effMassU, double effMassD, double effMassS, 
	double propagatorIntegralPrecision, scatteringProcess process, 
	bool largeAngleScatteringContribution, double crossSectionIntegralPrecision,
	double integralPrecision_dsdE, double integralPrecision_ds
)
{	
	IntegratedCrossSectionIntegrand aux(
		"integratedCrossSectionIntegral_dx", parametersNJL, T, 
        effChemPotU, effChemPotD, effChemPotS, 
        effMassU, effMassD, effMassS, 
        propagatorIntegralPrecision, process, 
        largeAngleScatteringContribution, crossSectionIntegralPrecision,
        integralPrecision_dsdE
	);

	double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMassU, effMassD, effMassS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effChemPotU, effChemPotD, effChemPotS, process, cP1, cP2);

	//for the integral over the Fermi-Dirac distributions, we use the same precision as for the integrations in the meson propagator
	double fermiDiracIntegral1 = fermiDiracIntegral(Nc, T, cP1, m1, propagatorIntegralPrecision);
    double fermiDiracIntegral2 = fermiDiracIntegral(Nc, T, cP2, m2, propagatorIntegralPrecision);

    double cte = 0.5*( ( Nc/(M_PI*M_PI) )*( Nc/(M_PI*M_PI) ) )/( fermiDiracIntegral1*fermiDiracIntegral2 );

    double sMin = pow(m1+m2,2);
    double sMax = pow(sqrt( pow(cutoff,2) + pow(m1,2) ) + sqrt( pow(cutoff,2) + pow(m2,2) ), 2);

    Integration1DimNewtonCotes trapezoidalSum(1.001*sMin, 0.999*sMax, 10, &aux, integratedCrossSectionIntegrand_ds, alternativeCompositeSimpson);
    double normalization_ds = trapezoidalSum.evaluate();
    aux.setNormalizationRiemannSum_ds(normalization_ds);

    int integrationWorkspace = 1000;
    Integration1DimGSLQAGS integratedCrossSectionIntegral_dx(1.0/( 1.0 + sMax ), 1.0/( 1.0 + sMin ), &aux, integratedCrossSectionCOVNormalizedIntegrand_dx, integralPrecision_ds, integralPrecision_ds, integrationWorkspace);
    double integratedCrossSection = normalization_ds*cte*(1.0/4.0)*integratedCrossSectionIntegral_dx.evaluate();

    return integratedCrossSection;
}


double integratedCrossSectionOGIntegrand_dp1dp2dtheta(double theta, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	double effCPD = aux.getDownQuarkEffectiveChemicalPotential();
	double effCPS = aux.getStrangeQuarkEffectiveChemicalPotential();
	double effMU = aux.getUpQuarkEffectiveMass();
	double effMD = aux.getDownQuarkEffectiveMass();
	double effMS = aux.getStrangeQuarkEffectiveMass();
	double p1 = aux.getMomentumParticle1();
	double p2 = aux.getMomentumParticle2();
	scatteringProcess process = aux.getProcess();

	//set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMU, effMD, effMS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effCPU, effCPD, effCPS, process, cP1, cP2);

	double E1 = sqrt( pow(p1,2) + pow(m1,2) );
	double E2 = sqrt( pow(p2,2) + pow(m2,2) );
	double fermiDist1 = fermiDistribution(T, E1 - cP1);
	double fermiDist2 = fermiDistribution(T, E2 - cP2);

	double s = pow(E1+E2, 2) - ( pow(p1,2) + pow(p2,2) + 2.0*p1*p2*cos(theta) );
	double vRel = momentumCM(s, m1, m2)*sqrt(s)/( E1*E2 );
	
	double crossSection = 1.0;
	crossSection = crossSectionProcess12To34(
		aux.getParametersNJL(), T, 
        effCPU, effCPD, effCPS, 
        effMU, effMD, effMS, 
        s, aux.getPropagatorIntegralPrecision(), process, 
        aux.getLargeAngleScatteringContribution(), aux.getCrossSectionIntegralPrecision()
	);

	double integratedCrossSection_dp1dp2dtheta = ( sin(theta)*pow(p1,2)*pow(p2,2) )*( fermiDist1*fermiDist2*vRel*crossSection );

    return integratedCrossSection_dp1dp2dtheta;
}


double integratedCrossSectionOGIntegrand_dp1dp2(double p2, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	scatteringProcess process = aux.getProcess();
	double integralPrecision = aux.getIntegratedCrossSectionIntegralPrecision_dXdYdZ();

	//set momentum particle 2
	aux.setMomentumParticle2(p2);

	double integratedCrossSection_dp1dp2 = 0.0;

	int integrationWorkspace = 1000;
    aux.setIntegralID("integratedCrossSectionOGIntegral_dp1dp2dtheta");
    Integration1DimGSLQAGS integratedCrossSectionOGIntegral_dp1dp2dtheta(0, M_PI, &aux, integratedCrossSectionOGIntegrand_dp1dp2dtheta, integralPrecision, integralPrecision, integrationWorkspace);
    integratedCrossSection_dp1dp2 = integratedCrossSection_dp1dp2 + integratedCrossSectionOGIntegral_dp1dp2dtheta.evaluate();

	cout << toString(process) << "\t" 
		 << "T = " << T << "\t" 
		 << "cPU = " << effCPU << "\t" 
		 << "p2 = " << p2 << "\t" <<  "integratedCrossSection_dp1dp2 = " << integratedCrossSection_dp1dp2 << "\n";

    return integratedCrossSection_dp1dp2;
}


double integratedCrossSectionOGIntegrand_dp1(double p1, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double cutoff = aux.getParametersNJL().getThreeMomentumCutoff();
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	scatteringProcess process = aux.getProcess();
	double integralPrecision = aux.getIntegratedCrossSectionIntegralPrecision_dXdY();

	//set momentum particle 1
	aux.setMomentumParticle1(p1);

	double integratedCrossSection_dp1 = 0.0;

	int integrationWorkspace = 1000;
    aux.setIntegralID("integratedCrossSectionOGIntegral_dp1dp2");
    Integration1DimGSLQAGS integratedCrossSectionOGIntegral_dp1dp2(0, cutoff, &aux, integratedCrossSectionOGIntegrand_dp1dp2, integralPrecision, integralPrecision, integrationWorkspace);
    integratedCrossSection_dp1 = integratedCrossSection_dp1 + integratedCrossSectionOGIntegral_dp1dp2.evaluate();

	cout << toString(process) << "\t" 
		 << "T = " << T << "\t" 
		 << "cPU = " << effCPU << "\t" 
		 << "p1 = " << p1 << "\t" <<  "integratedCrossSection_dp1 = " << integratedCrossSection_dp1 << "\n";

    return integratedCrossSection_dp1;
}


double integratedCrossSectionOGProcess12To34(
	SU3NJL3DCutoffParameters parametersNJL, double T, 
	double effChemPotU, double effChemPotD, double effChemPotS, 
	double effMassU, double effMassD, double effMassS, 
	double propagatorIntegralPrecision, scatteringProcess process, 
	bool largeAngleScatteringContribution, double crossSectionIntegralPrecision,
	double integralPrecision_dp1dp2dtheta, double integralPrecision_dp1dp2, double integralPrecision_dp1
)
{   
	IntegratedCrossSectionIntegrand aux(
		"integratedCrossSectionOGIntegral_dp1", parametersNJL, T, 
        effChemPotU, effChemPotD, effChemPotS, 
        effMassU, effMassD, effMassS, 
        propagatorIntegralPrecision, process, 
        largeAngleScatteringContribution, crossSectionIntegralPrecision,
        integralPrecision_dp1dp2dtheta, integralPrecision_dp1dp2
	);

	double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMassU, effMassD, effMassS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effChemPotU, effChemPotD, effChemPotS, process, cP1, cP2);

	//for the integral over the Fermi-Dirac distributions, we use the same precision as for the integrations in the meson propagator
	double fermiDiracIntegral1 = fermiDiracIntegral(Nc, T, cP1, m1, propagatorIntegralPrecision);
    double fermiDiracIntegral2 = fermiDiracIntegral(Nc, T, cP2, m2, propagatorIntegralPrecision);

    double cte = 0.5*( ( Nc/(M_PI*M_PI) )*( Nc/(M_PI*M_PI) ) )/( fermiDiracIntegral1*fermiDiracIntegral2 );

    int integrationWorkspace = 1000;
    Integration1DimGSLQAGS integratedCrossSectionOGIntegral_dp1(0.0, cutoff, &aux, integratedCrossSectionOGIntegrand_dp1, integralPrecision_dp1, integralPrecision_dp1, integrationWorkspace);
    double integratedCrossSection = cte*integratedCrossSectionOGIntegral_dp1.evaluate();

    return integratedCrossSection;
}


////////////////////////////////////////////////////////////////////////////////////////
//12->34 average transition rate W Klevansky approximation: integrating up to the cutoff
//temperature and chemical potential dependent probability (a la Klevansky)
// to produce a quark-quark or quark-antiquark pair in the medium 


double probabilityKlevanskyHeavisideSolutionPlus(double m1, double m2, double s, double E1)
{
	double sqrtArg = ( pow(E1,2) - pow(m1,2) )*( pow(m1,4) + pow(pow(m2,2)-s,2) - 2.0*pow(m1,2)*( pow(m2,2) + s ) );

	double aux = ( E1*( s - pow(m1,2) - pow(m2,2) ) + sqrt( sqrtArg ) )/( 2.0*pow(m1,2) );

	return aux;
}


double probabilityKlevanskyHeavisideSolutionMinus(double m1, double m2, double s, double E1)
{
	double sqrtArg = ( pow(E1,2) - pow(m1,2) )*( pow(m1,4) + pow(pow(m2,2)-s,2) - 2.0*pow(m1,2)*( pow(m2,2) + s ) );

	double aux = ( E1*( s - pow(m1,2) - pow(m2,2) ) - sqrt( sqrtArg ) )/( 2.0*pow(m1,2) );

	return aux;
}


double probabilityKlevanskyIntegrand_dx(double x, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	double effCPD = aux.getDownQuarkEffectiveChemicalPotential();
	double effCPS = aux.getStrangeQuarkEffectiveChemicalPotential();
	double effMU = aux.getUpQuarkEffectiveMass();
	double effMD = aux.getDownQuarkEffectiveMass();
	double effMS = aux.getStrangeQuarkEffectiveMass();
	double s = aux.getCenterOfMassEnergy();
	scatteringProcess process = aux.getProcess();

	//set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMU, effMD, effMS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effCPU, effCPD, effCPS, process, cP1, cP2);

	//change of variables
	double E1 = ( 1-x )/x;

    double fermiDist1 = fermiDistribution(T, E1 - cP1);
    double E2Minus = probabilityKlevanskyHeavisideSolutionMinus(m1, m2, s, E1);
    double E2Plus = probabilityKlevanskyHeavisideSolutionPlus(m1, m2, s, E1);
	double nonNormalizedKlevanskyIntegrand = fermiDist1*( log( 1.0 + exp(-(E2Minus-cP2)/T) ) - log( 1.0 + exp(-(E2Plus-cP2)/T) ) );

	nonNormalizedKlevanskyIntegrand = nonNormalizedKlevanskyIntegrand/pow(x,2);//needed because of the change of variables

    return nonNormalizedKlevanskyIntegrand;
}


double probabilityKlevansky(
	SU3NJL3DCutoffParameters parametersNJL, double T, 
	double effChemPotU, double effChemPotD, double effChemPotS, 
	double effMassU, double effMassD, double effMassS, 
	double s, scatteringProcess process,
	double integralPrecision_dE1
)
{   
	IntegratedCrossSectionIntegrand aux(
		"probabilityKlevanskyIntegral_dx", parametersNJL, T, 
        effChemPotU, effChemPotD, effChemPotS, 
    	effMassU, effMassD, effMassS, 
        0.0, process, 
        false, 0.0,
        0.0
	);

	//set center of mass energy
	aux.setCenterOfMassEnergy(s);

    //set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMassU, effMassD, effMassS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effChemPotU, effChemPotD, effChemPotS, process, cP1, cP2);

	//get the number of colours
	double Nc = aux.getParametersNJL().getNumberOfColours();

	int integrationWorkspace = 1000;
	aux.setIntegralID("probabilityKlevanskyIntegral_dx");
    Integration1DimGSLQAGS probabilityKlevanskyIntegral_dx(0.0, 1.0/( 1.0 + m1 ), &aux, probabilityKlevanskyIntegrand_dx, integralPrecision_dE1, integralPrecision_dE1, integrationWorkspace);
    double probabilityKlev = (2.0*Nc)*(2.0*Nc)*(T)*probabilityKlevanskyIntegral_dx.evaluate();

    return probabilityKlev;
}


double probabilityKlevansky(double s, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	double effCPD = aux.getDownQuarkEffectiveChemicalPotential();
	double effCPS = aux.getStrangeQuarkEffectiveChemicalPotential();
	double effMU = aux.getUpQuarkEffectiveMass();
	double effMD = aux.getDownQuarkEffectiveMass();
	double effMS = aux.getStrangeQuarkEffectiveMass();
	scatteringProcess process = aux.getProcess();
	double integralPrecision_dsdx = aux.getIntegratedCrossSectionIntegralPrecision_dXdY();

	double probabilityKlev = probabilityKlevansky(
		aux.getParametersNJL(), T, 
		effCPU, effCPD, effCPS, 
		effMU, effMD, effMS, 
		s, process,
		integralPrecision_dsdx
	);

    return probabilityKlev;
}


double integratedCrossSectionKlevanskyIntegrand_ds(double s, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	double effCPD = aux.getDownQuarkEffectiveChemicalPotential();
	double effCPS = aux.getStrangeQuarkEffectiveChemicalPotential();
	double effMU = aux.getUpQuarkEffectiveMass();
	double effMD = aux.getDownQuarkEffectiveMass();
	double effMS = aux.getStrangeQuarkEffectiveMass();
	scatteringProcess process = aux.getProcess();

	//set center of mass energy
	aux.setCenterOfMassEnergy(s);

    //set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMU, effMD, effMS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effCPU, effCPD, effCPS, process, cP1, cP2);

	double probability = probabilityKlevansky(s, parameters);

	double crossSection = 1.0;
	crossSection = crossSectionProcess12To34(
		aux.getParametersNJL(), T, 
        effCPU, effCPD, effCPS, 
        effMU, effMD, effMS, 
        s, aux.getPropagatorIntegralPrecision(), process, 
        aux.getLargeAngleScatteringContribution(), aux.getCrossSectionIntegralPrecision()
	);

	double E1 = energyOutParticle1(s, m1, m2);
	double E2 = energyOutParticle2(s, m1, m2);
	double pCM = momentumCM(s, m1, m2);
	double integratedCrossSection = sqrt( pow(E1*E2 + pCM*pCM, 2) - pow(m1*m2, 2) )*crossSection*probability;

	cout << toString(process) << "\t" 
		 << "T = " << T << "\t" 
		 << "cPU = " << effCPU << "\t" 
		 << "s = " << s << "\t" <<  "integratedCrossSectionKlevansky = " << integratedCrossSection << "\n";

    return integratedCrossSection;
}


double integratedCrossSectionKlevanskyNormalizedIntegrand_dx(double x, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double normalization = aux.getNormalizationRiemannSum_ds();

    //change of variables
    double s = ( 1-x )/x;

    double integratedCrossSection_dx = integratedCrossSectionKlevanskyIntegrand_ds(s, parameters);

    integratedCrossSection_dx = integratedCrossSection_dx/pow(x,2);

    integratedCrossSection_dx = integratedCrossSection_dx/normalization;

    return integratedCrossSection_dx;
}


double integratedCrossSectionProcess12To34Klevansky(
	SU3NJL3DCutoffParameters parametersNJL, double T, 
	double effChemPotU, double effChemPotD, double effChemPotS, 
	double effMassU, double effMassD, double effMassS, 
	double propagatorIntegralPrecision, scatteringProcess process, 
	bool largeAngleScatteringContribution, double crossSectionIntegralPrecision,
	double integralPrecision_dsdE, double integralPrecision_ds
)
{	
	IntegratedCrossSectionIntegrand aux(
		"integratedCrossSectionKlevasnkyIntegral_dx", parametersNJL, T, 
        effChemPotU, effChemPotD, effChemPotS, 
    	effMassU, effMassD, effMassS, 
        propagatorIntegralPrecision, process, 
        largeAngleScatteringContribution, crossSectionIntegralPrecision,
        integralPrecision_dsdE
	);

	double cutoff = parametersNJL.getThreeMomentumCutoff();
	double Nc = parametersNJL.getNumberOfColours();

	//set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMassU, effMassD, effMassS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effChemPotU, effChemPotD, effChemPotS, process, cP1, cP2);

	//for the integral over the Fermi-Dirac distributions, we use the same precision as for the integrations in the meson propagator
	double fermiDiracIntegral1 = fermiDiracIntegral(Nc, T, cP1, m1, propagatorIntegralPrecision);
    double fermiDiracIntegral2 = fermiDiracIntegral(Nc, T, cP2, m2, propagatorIntegralPrecision);

    double sMin = centerOfMassEnergyThreshold(m1, m2, m3, m4);
    double sMax = sMaximumKlevansky(cutoff, effMassU, effMassD, effMassS);

    Integration1DimNewtonCotes trapezoidalSum(1.001*sMin, 0.999*sMax, 10, &aux, integratedCrossSectionKlevanskyIntegrand_ds, alternativeCompositeSimpson);
    double normalization_ds = trapezoidalSum.evaluate();
    aux.setNormalizationRiemannSum_ds(normalization_ds);

    int integrationWorkspace = 1000;
    Integration1DimGSLQAGS integratedCrossSectionKlevasnkyIntegral_dx(1.0/( 1.0 + sMax ), 1.0/( 1.0 + sMin ), &aux, integratedCrossSectionKlevanskyNormalizedIntegrand_dx, integralPrecision_ds, integralPrecision_ds, integrationWorkspace);
    double integratedCrossSection = integratedCrossSectionKlevasnkyIntegral_dx.evaluate();
    integratedCrossSection = normalization_ds*integratedCrossSection/( fermiDiracIntegral1*fermiDiracIntegral2*16*(M_PI*M_PI*M_PI*M_PI) );

    return integratedCrossSection;
}


////////////////////////////////////////////////////////////////////////////////////////
//q1q2->q3q4 integrated cross section integrating up to the cutoff using the Zhuang approximation


//temperature and chemical potential dependent constant that normalizes the probability to produce a quark-quark or quark-antiquark pair in the medium 
double nonNormalizedProbabilityZhuangIntegrand_ds(double s, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	double effCPD = aux.getDownQuarkEffectiveChemicalPotential();
	double effCPS = aux.getStrangeQuarkEffectiveChemicalPotential();
	double effMU = aux.getUpQuarkEffectiveMass();
	double effMD = aux.getDownQuarkEffectiveMass();
	double effMS = aux.getStrangeQuarkEffectiveMass();
	scatteringProcess process = aux.getProcess();

	//set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMU, effMD, effMS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effCPU, effCPD, effCPS, process, cP1, cP2);

	double pCM = momentumCM(s, m1, m2);
	double E1 = energyOutParticle1(s, m1, m2);
	double E2 = energyOutParticle2(s, m1, m2);
	double vRel = momentumCM(s, m1, m2)*sqrt(s)/( E1*E2 );
	double fermiDist1 = fermiDistribution(T, E1 - cP1);
	double fermiDist2 = fermiDistribution(T, E2 - cP2);

	double normalization = E1*E2*pCM*vRel*fermiDist1*fermiDist2;

    return normalization;
}


double probabilityNormalizationInverseZhuang(
	SU3NJL3DCutoffParameters parametersNJL, double T, 
	double effChemPotU, double effChemPotD, double effChemPotS, 
	double effMassU, double effMassD, double effMassS, 
	scatteringProcess process,
	double integralPrecision_ds
)
{
	IntegratedCrossSectionIntegrand aux(
		"nonNormalizedProbabilityZhuangIntegral_ds", parametersNJL, T, 
        effChemPotU, effChemPotD, effChemPotS, 
        effMassU, effMassD, effMassS, 
        0.0, process, 
        false, 0.0,
    	0.0
	);

    //set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMassU, effMassD, effMassS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effChemPotU, effChemPotD, effChemPotS, process, cP1, cP2);

    double sMin = centerOfMassEnergyThreshold(m1, m2, m3, m4);

    int integrationWorkspace = 1000;
    Integration1DimGSLQAGIU nonNormalizedProbabilityZhuangIntegral_ds(sMin, &aux, nonNormalizedProbabilityZhuangIntegrand_ds, integralPrecision_ds, integralPrecision_ds, integrationWorkspace);
    double normalization = nonNormalizedProbabilityZhuangIntegral_ds.evaluate();

    return 1.0/normalization;
}


double integratedCrossSectionZhuangIntegrand_ds(double s, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double T = aux.getTemperature();
	double effCPU = aux.getUpQuarkEffectiveChemicalPotential();
	double effCPD = aux.getDownQuarkEffectiveChemicalPotential();
	double effCPS = aux.getStrangeQuarkEffectiveChemicalPotential();
	double effMU = aux.getUpQuarkEffectiveMass();
	double effMD = aux.getDownQuarkEffectiveMass();
	double effMS = aux.getStrangeQuarkEffectiveMass();
	scatteringProcess process = aux.getProcess();

	//set center of mass energy
	aux.setCenterOfMassEnergy(s);

    //set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMU, effMD, effMS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effCPU, effCPD, effCPS, process, cP1, cP2);

	double nonNormalizedProbability = nonNormalizedProbabilityZhuangIntegrand_ds(s, parameters);

	double crossSection = 1.0;
	crossSection = crossSectionProcess12To34(
		aux.getParametersNJL(), T, 
        effCPU, effCPD, effCPS, 
        effMU, effMD, effMS, 
        s, aux.getPropagatorIntegralPrecision(), process, 
        aux.getLargeAngleScatteringContribution(), aux.getCrossSectionIntegralPrecision()
	);

	double integratedCrossSection = crossSection*nonNormalizedProbability;

	cout << toString(process) << "\t" 
		 << "T = " << T << "\t" 
		 << "cPU = " << effCPU << "\t" 
		 << "s = " << s << "\t" <<  "integratedCrossSectionZhuang = " << integratedCrossSection << "\n";

    return integratedCrossSection;
}


double integratedCrossSectionZhuangNormalizedIntegrand_dx(double x, void *parameters)
{   
	IntegratedCrossSectionIntegrand aux(parameters);
	double normalization = aux.getNormalizationRiemannSum_ds();

    //change of variables
    double s = ( 1-x )/x;

    double integratedCrossSection_dx = integratedCrossSectionZhuangIntegrand_ds(s, parameters);

    integratedCrossSection_dx = integratedCrossSection_dx/pow(x,2);

    integratedCrossSection_dx = integratedCrossSection_dx/normalization;

    return integratedCrossSection_dx;
}


double integratedCrossSectionProcess12To34Zhuang(
	SU3NJL3DCutoffParameters parametersNJL, double T, 
	double effChemPotU, double effChemPotD, double effChemPotS, 
	double effMassU, double effMassD, double effMassS, 
	double propagatorIntegralPrecision, scatteringProcess process, 
	bool largeAngleScatteringContribution, double crossSectionIntegralPrecision,
	double integralPrecision_ds
)
{	
	IntegratedCrossSectionIntegrand aux(
		"integratedCrossSectionZhuangNormalizedIntegral_dx", parametersNJL, T, 
        effChemPotU, effChemPotD, effChemPotS, 
        effMassU, effMassD, effMassS, 
        propagatorIntegralPrecision, process, 
        largeAngleScatteringContribution, crossSectionIntegralPrecision,
    	0.0
	);

	double cutoff = parametersNJL.getThreeMomentumCutoff();

	//set incoming and outgoing masses
	double m1, m2, m3, m4;
	inOutMassesGivenScatteringProcess(effMassU, effMassD, effMassS, process, m1, m2, m3, m4);

	//set incoming chemical potentials
	double cP1, cP2;
	inChemicalPotentialsGivenScatteringProcess(effChemPotU, effChemPotD, effChemPotS, process, cP1, cP2);

    double sMin = centerOfMassEnergyThreshold(m1, m2, m3, m4);
    double sMax = sMaximumKlevansky(cutoff, effMassU, effMassD, effMassS);

    double cte = probabilityNormalizationInverseZhuang(
		parametersNJL, T,
		effChemPotU, effChemPotD, effChemPotS, 
		effMassU, effMassD, effMassS, 
		process,
		integralPrecision_ds
	);

    Integration1DimNewtonCotes trapezoidalSum(1.001*sMin, 0.999*sMax, 10, &aux, integratedCrossSectionZhuangIntegrand_ds, alternativeCompositeSimpson);
    double normalization_ds = trapezoidalSum.evaluate();
    aux.setNormalizationRiemannSum_ds(normalization_ds);

	int integrationWorkspace = 1000;
    Integration1DimGSLQAGS integratedCrossSectionZhuangNormalizedIntegral_dx(1.0/( 1.0 + sMax ), 1.0/( 1.0 + sMin ), &aux, integratedCrossSectionZhuangNormalizedIntegrand_dx, integralPrecision_ds, integralPrecision_ds, integrationWorkspace);
    double integratedCrossSection = normalization_ds*cte*integratedCrossSectionZhuangNormalizedIntegral_dx.evaluate();

    return integratedCrossSection;
}


string toString(IntegratedCrossSectionApproximationMethod method) 
{
    // Check if the method exists in the map using count
    if (IntegratedCrossSectionApproximationMethodMap.count(method))
    {
        return IntegratedCrossSectionApproximationMethodMap.at(method);
    } 
    else 
    {
        cout << "Error: IntegratedCrossSectionApproximationMethod not found in map! Returning UNKNOWN." << endl;
        return "UNKNOWN";
    }
}

IntegratedCrossSectionApproximationMethod stringToIntegratedCrossSectionApproximationMethod(const std::string& methodString) 
{
    // Iterate over the map with explicit type
    for (map<IntegratedCrossSectionApproximationMethod, string>::const_iterator it = IntegratedCrossSectionApproximationMethodMap.begin(); it != IntegratedCrossSectionApproximationMethodMap.end(); ++it) 
    {
        if (it->second == methodString) 
        {
            return it->first;
        }
    }

    std::cout << "Invalid IntegratedCrossSectionApproximationMethod string: " << methodString << ". Aborting!" << std::endl;
    abort();
}

bool isValidIntegratedCrossSectionApproximationMethod(const string& methodString)
{
    bool isIntegratedCrossSectionApproximationMethod = false;
    // Iterate over the map with explicit type
    for (map<IntegratedCrossSectionApproximationMethod, string>::const_iterator it = IntegratedCrossSectionApproximationMethodMap.begin(); it != IntegratedCrossSectionApproximationMethodMap.end(); ++it) 
    {
        if (it->second == methodString) 
        {
            isIntegratedCrossSectionApproximationMethod = true;
            break;
        }
    }

    if( isIntegratedCrossSectionApproximationMethod==false )
    {
        cout << "The value " + methodString + " is not a IntegratedCrossSectionApproximationMethod!\n";
    }

    return isIntegratedCrossSectionApproximationMethod;
}

vector<SU3NJL3DCutoffIntegratedCrossSection> evaluateIntegratedCrossSectionAlongTrajectory(
	vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTemperatureSolution,
	scatteringProcess process, 
	double propagatorIntegralPrecision,
	bool largeAngleScatteringContribution, 
	double crossSectionIntegralPrecision,
	double integratedCrossSectionIntegralPrecision_dXdY, 
	double integratedCrossSectionIntegralPrecision_dX, 
	IntegratedCrossSectionApproximationMethod approximationMethod,
	int numberOfThreads
)
{	
	//get parameters from provided set of solutions
	SU3NJL3DCutoffParameters parameters = finiteTemperatureSolution[0].getParametersNJL();

	//vector to save integrated cross sections
    int size = int(finiteTemperatureSolution.size());
    vector<SU3NJL3DCutoffIntegratedCrossSection> integratedCrossSectionFiniteTemperature(size);

    cout << "Number of threads being used: " << numberOfThreads << "\n";
    #pragma omp parallel for schedule(dynamic) num_threads( numberOfThreads )
    for (int i = 0; i < size; ++i)
    {
        double T = finiteTemperatureSolution[i].getTemperature();
        double effChemPotU = finiteTemperatureSolution[i].getUpQuarkChemicalPotential();
        double effChemPotD = finiteTemperatureSolution[i].getDownQuarkChemicalPotential();
        double effChemPotS = finiteTemperatureSolution[i].getStrangeQuarkChemicalPotential();
        double effMassU = finiteTemperatureSolution[i].getUpQuarkEffectiveMass();
        double effMassD = finiteTemperatureSolution[i].getDownQuarkEffectiveMass();
        double effMassS = finiteTemperatureSolution[i].getStrangeQuarkEffectiveMass();

        SU3NJL3DCutoffIntegratedCrossSection intCrossSection(
			parameters, T, 
            effChemPotU, effChemPotD, effChemPotS, 
            effMassU, effMassD, effMassS, 
            propagatorIntegralPrecision, process,
            largeAngleScatteringContribution, crossSectionIntegralPrecision,
            integratedCrossSectionIntegralPrecision_dXdY, integratedCrossSectionIntegralPrecision_dX,
            approximationMethod
		);
        intCrossSection.setIntegratedCrossSection();
        intCrossSection.setQuarkNumbers();

        integratedCrossSectionFiniteTemperature[i] = intCrossSection;
    }

    return integratedCrossSectionFiniteTemperature;
}


void writeIntegratedCrossSectionToFile(vector<SU3NJL3DCutoffIntegratedCrossSection> integratedCrossSection, string fileName)
{
    std::ofstream fileTest;
    fileTest.open(fileName, std::ofstream::out | std::ios::trunc);
    fileTest.precision(15);
    fileTest.width(25);   fileTest << "T[GeV]";           //temperature
    fileTest.width(25);   fileTest << "nU[GeV^3]";       //u quark number
    fileTest.width(25);   fileTest << "nD[GeV^3]";       //d quark number
    fileTest.width(25);   fileTest << "nS[GeV^3]";       //s quark number
    fileTest.width(25);   fileTest << "nUBar[GeV^3]";    //u anti-quark number
    fileTest.width(25);   fileTest << "nDBar[GeV^3]";    //d anti-quark number
    fileTest.width(25);   fileTest << "nSBbar[GeV^3]";    //s anti-quark number
    fileTest.width(25);   fileTest << "W[Gev^-2]";        //integrated cross section
    fileTest.width(25);   fileTest << "effMU[GeV]";          //u-quark effective mass
    fileTest.width(25);   fileTest << "effMD[GeV]";          //d-quark effective mass
    fileTest.width(25);   fileTest << "effMS[GeV]";          //s-quark effective mass
    fileTest.width(25);   fileTest << "effCPU[GeV]";          //u quark effective chemical potential
    fileTest.width(25);   fileTest << "effCPD[GeV]";          //d quark effective chemical potential
    fileTest.width(25);   fileTest << "effCPS[GeV]";          //s quark effective chemical potential
    fileTest << std::endl;

    for (int i = 0; i < int(integratedCrossSection.size()); ++i)
    {
        fileTest.width(25);   fileTest << integratedCrossSection[i].getTemperature(); //temperature
        fileTest.width(25);   fileTest << integratedCrossSection[i].getUpQuarkNumber(); //u-quark number
        fileTest.width(25);   fileTest << integratedCrossSection[i].getDownQuarkNumber(); //d-quark number
        fileTest.width(25);   fileTest << integratedCrossSection[i].getStrangeQuarkNumber(); //s-quark number
        fileTest.width(25);   fileTest << integratedCrossSection[i].getUpAntiquarkNumber(); //u-antiquark number
        fileTest.width(25);   fileTest << integratedCrossSection[i].getDownAntiquarkNumber(); //d-antiquark number
        fileTest.width(25);   fileTest << integratedCrossSection[i].getStrangeAntiquarkNumber(); //s-antiquark number
        fileTest.width(25);   fileTest << integratedCrossSection[i].getIntegratedCrossSection(); //integrated cross section
        fileTest.width(25);   fileTest << integratedCrossSection[i].getUpQuarkEffectiveMass(); //u-quark effective mass
        fileTest.width(25);   fileTest << integratedCrossSection[i].getDownQuarkEffectiveMass(); //d-quark effective mass
        fileTest.width(25);   fileTest << integratedCrossSection[i].getStrangeQuarkEffectiveMass(); //s-quark effective mass
        fileTest.width(25);   fileTest << integratedCrossSection[i].getUpQuarkEffectiveChemicalPotential(); //u-quark chemical potential
        fileTest.width(25);   fileTest << integratedCrossSection[i].getDownQuarkEffectiveChemicalPotential(); //u-quark chemical potential
        fileTest.width(25);   fileTest << integratedCrossSection[i].getStrangeQuarkEffectiveChemicalPotential(); //u-quark chemical potential
        fileTest << std::endl;
    }

    fileTest.close();
}


void evaluateIntegratedCrossSectionAlongFixedChemicalPotentialTrajectory(
	vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTemperatureSolution,
	scatteringProcess process, 
	double propagatorIntegralPrecision,
	bool largeAngleScatteringContribution, 
	double crossSectionIntegralPrecision,
	double integratedCrossSectionIntegralPrecision_dXdY, 
	double integratedCrossSectionIntegralPrecision_dX, 
	IntegratedCrossSectionApproximationMethod approximationMethod,
	int numberOfThreads
)
{
	//calculate integrated cross sections for the path at finite temperature
    vector<SU3NJL3DCutoffIntegratedCrossSection> integratedCrossSectionFiniteTemperature =
    evaluateIntegratedCrossSectionAlongTrajectory(
		finiteTemperatureSolution,
        process, 
        propagatorIntegralPrecision,
        largeAngleScatteringContribution, 
        crossSectionIntegralPrecision,
        integratedCrossSectionIntegralPrecision_dXdY, 
        integratedCrossSectionIntegralPrecision_dX, 
        approximationMethod,
		numberOfThreads
	);


    string fileName = "IntegratedCrossSection";
    fileName = fileName + "_" + integratedCrossSectionFiniteTemperature[0].getParametersNJL().getParameterSetName();
    fileName = fileName + "_" + toString(integratedCrossSectionFiniteTemperature[0].getProcess());
    if ( integratedCrossSectionFiniteTemperature[0].getLargeAngleScatteringContribution() ){ fileName = fileName + "_LargeAngleScat"; }
    fileName = fileName + "_" + toString(approximationMethod);
    fileName = fileName + "_TMin" + to_string(integratedCrossSectionFiniteTemperature[0].getTemperature());
    fileName = fileName + "_TMax" + to_string(integratedCrossSectionFiniteTemperature[integratedCrossSectionFiniteTemperature.size()-1].getTemperature());
    fileName = fileName + "_CPU" + to_string(integratedCrossSectionFiniteTemperature[0].getUpQuarkEffectiveChemicalPotential());
    //std::replace( fileName.begin(), fileName.end(), '.', ','); 
    fileName =  fileName +".dat";

    writeIntegratedCrossSectionToFile(integratedCrossSectionFiniteTemperature, fileName);
}


void evaluateIsospinSymmetricIntegratedCrossSectionsAlongFixedChemicalPotentialTrajectory(
	vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTemperatureSolution, 
	double propagatorIntegralPrecision, 
	bool largeAngleScatteringContribution, 
	double crossSectionIntegralPrecision, 
	double integratedCrossSectionIntegralPrecision_dXdY, 
	double integratedCrossSectionIntegralPrecision_dX, 
	IntegratedCrossSectionApproximationMethod approximationMethod,
	int numberOfThreads
)
{
	bool runFiniteDensityProcesses = false;
	for (int i = 0; i < int(finiteTemperatureSolution.size()); ++i)
	{
		double effCPU = finiteTemperatureSolution[i].getUpQuarkChemicalPotential();
		double effCPD = finiteTemperatureSolution[i].getUpQuarkChemicalPotential();
		double effCPS = finiteTemperatureSolution[i].getUpQuarkChemicalPotential();

		if ( effCPU>0 || effCPD>0 || effCPS>0 )
		{
			runFiniteDensityProcesses = true;
		}
	}

	vector<scatteringProcess> process;
	if ( runFiniteDensityProcesses )
	{
		process = { 
			UUUU, UDUD, USUS, SSSS,
			UUBarUUBar, UUBarDDBar, UUBarSSBar, UDBarUDBar, USBarUSBar, SUBarSUBar, SSBarUUBar, SSBarSSBar,
			UBarUBarUBarUBar, UBarDBarUBarDBar, UBarSBarUBarSBar, SBarSBarSBarSBar 
		};
	}
	else
	{
		process = { 
			UUUU, UDUD, USUS, SSSS,
			UUBarUUBar, UUBarDDBar, UUBarSSBar, UDBarUDBar, USBarUSBar, SSBarUUBar, SSBarSSBar 
		};
	}

	for (int i = 0; i < int(process.size()); ++i)
	{
		evaluateIntegratedCrossSectionAlongFixedChemicalPotentialTrajectory(
			finiteTemperatureSolution, 
			process[i], 
			propagatorIntegralPrecision, 
			largeAngleScatteringContribution, 
			crossSectionIntegralPrecision, 
			integratedCrossSectionIntegralPrecision_dXdY, 
			integratedCrossSectionIntegralPrecision_dX, 
			approximationMethod,
			numberOfThreads
		);
	}
}


void evaluateIntegratedCrossSectionAlongFixedTemperatureTrajectory(
	vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTemperatureSolution,
	scatteringProcess process, 
	double propagatorIntegralPrecision,
	bool largeAngleScatteringContribution, 
	double crossSectionIntegralPrecision,
	double integratedCrossSectionIntegralPrecision_dXdY, 
	double integratedCrossSectionIntegralPrecision_dX, 
	IntegratedCrossSectionApproximationMethod approximationMethod,
	int numberOfThreads
)
{
	//calculate integrated cross sections for the path at finite temperature
    vector<SU3NJL3DCutoffIntegratedCrossSection> integratedCrossSectionFiniteTemperature =
    evaluateIntegratedCrossSectionAlongTrajectory(
		finiteTemperatureSolution,
        process, 
        propagatorIntegralPrecision,
        largeAngleScatteringContribution, 
        crossSectionIntegralPrecision,
        integratedCrossSectionIntegralPrecision_dXdY, 
        integratedCrossSectionIntegralPrecision_dX, 
        approximationMethod,
		numberOfThreads
	);


    string fileName = "IntegratedCrossSection";
    fileName = fileName + "_" + integratedCrossSectionFiniteTemperature[0].getParametersNJL().getParameterSetName();
    fileName = fileName + "_" + toString(integratedCrossSectionFiniteTemperature[0].getProcess());
    if ( integratedCrossSectionFiniteTemperature[0].getLargeAngleScatteringContribution() ){ fileName = fileName + "_LargeAngleScat"; }
    fileName = fileName + "_" + toString(approximationMethod);
    fileName = fileName + "_CPQMin" + to_string(integratedCrossSectionFiniteTemperature[0].getUpQuarkEffectiveChemicalPotential());
    fileName = fileName + "_CPQMax" + to_string(integratedCrossSectionFiniteTemperature[integratedCrossSectionFiniteTemperature.size()-1].getUpQuarkEffectiveChemicalPotential());
    fileName = fileName + "_T" + to_string(integratedCrossSectionFiniteTemperature[0].getTemperature());
    //std::replace( fileName.begin(), fileName.end(), '.', ','); 
    fileName =  fileName +".dat";

    writeIntegratedCrossSectionToFile(integratedCrossSectionFiniteTemperature, fileName);
}


void evaluateIsospinSymmetricIntegratedCrossSectionsAlongFixedTemperatureTrajectory(
	vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTemperatureSolution,
	double propagatorIntegralPrecision,
	bool largeAngleScatteringContribution, 
	double crossSectionIntegralPrecision,
	double integratedCrossSectionIntegralPrecision_dXdY, 
	double integratedCrossSectionIntegralPrecision_dX, 
	IntegratedCrossSectionApproximationMethod approximationMethod,
	int numberOfThreads
)
{
	vector<scatteringProcess> process;
	process = { 
		UUUU, UDUD, USUS, SSSS,
		UUBarUUBar, UUBarDDBar, UUBarSSBar, UDBarUDBar, USBarUSBar, SUBarSUBar, SSBarUUBar, SSBarSSBar,
		UBarUBarUBarUBar, UBarDBarUBarDBar, UBarSBarUBarSBar, SBarSBarSBarSBar 
	};

	for (int i = 0; i < int(process.size()); ++i)
	{
		evaluateIntegratedCrossSectionAlongFixedTemperatureTrajectory(
			finiteTemperatureSolution,
        	process[i], 
            propagatorIntegralPrecision,
            largeAngleScatteringContribution, 
            crossSectionIntegralPrecision,
            integratedCrossSectionIntegralPrecision_dXdY, 
            integratedCrossSectionIntegralPrecision_dX, 
            approximationMethod,
			numberOfThreads
		);
	}
}


void evaluateIsospinSymmetricIntegratedCrossSectionsWithZeroChemicalPotential(
	SU3NJL3DCutoffParameters& parameters,                                    
    double precisionVacuum,                                    
    MultiRootFindingMethod methodVacuum,                                    
    double lightQuarkMassGuess, 
    double strangeQuarkMassGuess,
	double precisionVacToFinTemp,
	MultiRootFindingMethod methodVacToFinTemp, 
	double minimumTemperature, 
	double maximumTemperature, 
	int numberOfPointsFromVacToMinTemp, 
	int numberOfPointsFromMinToMaxTemp, 
	bool largeAngleScatteringContribution, 
	IntegratedCrossSectionApproximationMethod approximationMethod, 
	double propagatorIntegralPrecision, 
	double crossSectionIntegralPrecision, 
	double integratedCrossSectionIntegralPrecision_dXdY, 
	double integratedCrossSectionIntegralPrecision_dX,
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

    //solve model at zero chemical potential up to some finite temperature
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTSolution = 
    solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, minimumTemperature, numberOfPointsFromVacToMinTemp, precisionVacToFinTemp, methodVacToFinTemp);
    for (int i = 0; i < int(finiteTSolution.size()); ++i)
    {
        cout << finiteTSolution[i].getTemperature() << "\t" 
        	 << finiteTSolution[i].getUpQuarkEffectiveMass() << "\t"
        	 << finiteTSolution[i].getDownQuarkEffectiveMass() << "\t"
        	 << finiteTSolution[i].getStrangeQuarkEffectiveMass() << "\n";
    }

    finiteTSolution = 
    solveFromLowToHighTemperatureAtZeroChemicalPotential(finiteTSolution[finiteTSolution.size()-1], maximumTemperature, numberOfPointsFromMinToMaxTemp, precisionVacToFinTemp, methodVacToFinTemp);
    for (int i = 0; i < int(finiteTSolution.size()); ++i)
    {
        cout << finiteTSolution[i].getTemperature() << "\t" 
        	 << finiteTSolution[i].getUpQuarkEffectiveMass() << "\t"
        	 << finiteTSolution[i].getDownQuarkEffectiveMass() << "\t"
        	 << finiteTSolution[i].getStrangeQuarkEffectiveMass() << "\n";
    }

    evaluateIsospinSymmetricIntegratedCrossSectionsAlongFixedChemicalPotentialTrajectory(
		finiteTSolution, 
		propagatorIntegralPrecision, 
		largeAngleScatteringContribution, 
		crossSectionIntegralPrecision, 
		integratedCrossSectionIntegralPrecision_dXdY, 
		integratedCrossSectionIntegralPrecision_dX, 
		approximationMethod,
		numberOfThreads
	);
}


void evaluateIntegratedCrossSectionsWithFixedTemperature(
	SU3NJL3DCutoffVacuum vacuum,
	double gapPrecision,
	double fixedTemperature, 
	int numberOfPointsFromVacToMinTemp, 
	int numberOfPointsFromMinTempToMinChemPot, 
	int numberOfPointsFromMinToMaxChemPot, 
	double minChemPot,
	double maxChemPot,
	bool largeAngleScatteringContribution, 
	IntegratedCrossSectionApproximationMethod approximationMethod,
	double propagatorIntegralPrecision,
	double crossSectionIntegralPrecision,
	double integratedCrossSectionIntegralPrecision_dXdY,
	double integratedCrossSectionIntegralPrecision_dX,
	int numberOfThreads
)
{
    //solve model at zero chemical potential up to some finite temperature
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTSolution = 
    solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, fixedTemperature, numberOfPointsFromVacToMinTemp, gapPrecision, HYBRIDS);

    for (int i = 0; i < int(finiteTSolution.size()); ++i)
    {
        cout << finiteTSolution[i].getTemperature() << "\t"
             << finiteTSolution[i].getUpQuarkChemicalPotential() << "\t"
             << finiteTSolution[i].getDownQuarkChemicalPotential() << "\t"
             << finiteTSolution[i].getStrangeQuarkChemicalPotential() << "\t"
             << finiteTSolution[i].getUpQuarkEffectiveMass() << "\t"
             << finiteTSolution[i].getDownQuarkEffectiveMass() << "\t"
             << finiteTSolution[i].getStrangeQuarkEffectiveMass() << "\n";
    }


    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteChemPotSolution;
    if ( minChemPot>0 )
    {
    	finiteChemPotSolution = solveFromFiniteTemperatureToFiniteChemicalPotential(finiteTSolution[finiteTSolution.size()-1], minChemPot, numberOfPointsFromMinTempToMinChemPot, gapPrecision, HYBRIDS);
        for (int i = 0; i < int(finiteChemPotSolution.size()); ++i)
	    {   
	        cout << finiteChemPotSolution[i].getTemperature() << "\t"
	             << finiteChemPotSolution[i].getUpQuarkChemicalPotential() << "\t"
	             << finiteChemPotSolution[i].getDownQuarkChemicalPotential() << "\t"
	             << finiteChemPotSolution[i].getStrangeQuarkChemicalPotential() << "\t"
	             << finiteChemPotSolution[i].getUpQuarkEffectiveMass() << "\t"
	             << finiteChemPotSolution[i].getDownQuarkEffectiveMass() << "\t"
	             << finiteChemPotSolution[i].getStrangeQuarkEffectiveMass() << "\n";
	    }

    	finiteChemPotSolution = solveFromFiniteTemperatureToFiniteChemicalPotential(finiteChemPotSolution[finiteChemPotSolution.size()-1], maxChemPot, numberOfPointsFromMinToMaxChemPot, gapPrecision, HYBRIDS);
	    for (int i = 0; i < int(finiteChemPotSolution.size()); ++i)
	    {   
	        cout << finiteChemPotSolution[i].getTemperature() << "\t"
	             << finiteChemPotSolution[i].getUpQuarkChemicalPotential() << "\t"
	             << finiteChemPotSolution[i].getDownQuarkChemicalPotential() << "\t"
	             << finiteChemPotSolution[i].getStrangeQuarkChemicalPotential() << "\t"
	             << finiteChemPotSolution[i].getUpQuarkEffectiveMass() << "\t"
	             << finiteChemPotSolution[i].getDownQuarkEffectiveMass() << "\t"
	             << finiteChemPotSolution[i].getStrangeQuarkEffectiveMass() << "\n";
	    }
    }
    else
    {
    	finiteChemPotSolution = solveFromFiniteTemperatureToFiniteChemicalPotential(finiteTSolution[finiteTSolution.size()-1], maxChemPot, numberOfPointsFromMinTempToMinChemPot, gapPrecision, HYBRIDS);
	    for (int i = 0; i < int(finiteChemPotSolution.size()); ++i)
	    {   
	        cout << finiteChemPotSolution[i].getTemperature() << "\t"
	             << finiteChemPotSolution[i].getUpQuarkChemicalPotential() << "\t"
	             << finiteChemPotSolution[i].getDownQuarkChemicalPotential() << "\t"
	             << finiteChemPotSolution[i].getStrangeQuarkChemicalPotential() << "\t"
	             << finiteChemPotSolution[i].getUpQuarkEffectiveMass() << "\t"
	             << finiteChemPotSolution[i].getDownQuarkEffectiveMass() << "\t"
	             << finiteChemPotSolution[i].getStrangeQuarkEffectiveMass() << "\n";
	    }
	}

    evaluateIsospinSymmetricIntegratedCrossSectionsAlongFixedTemperatureTrajectory(
		finiteChemPotSolution, 
        propagatorIntegralPrecision,
        largeAngleScatteringContribution, 
        crossSectionIntegralPrecision,
        integratedCrossSectionIntegralPrecision_dXdY, 
    	integratedCrossSectionIntegralPrecision_dX, 
        approximationMethod,
		numberOfThreads
	);

}


void evaluateIntegratedCrossSectionsWithFixedChemicalPotential(
	SU3NJL3DCutoffVacuum vacuum,
	double gapPrecision,
	double chemPot,
	double minimumTemperature, 
	double maximumTemperature, 
	int numberOfPointsFromVacToMinTemp, 
	int numberOfPointsMinTempToChemPot, 
	int numberOfPointsFromMinToMaxTemp, 
	bool largeAngleScatteringContribution, 
	IntegratedCrossSectionApproximationMethod approximationMethod,
	double propagatorIntegralPrecision,
	double crossSectionIntegralPrecision,
	double integratedCrossSectionIntegralPrecision_dXdY,
	double integratedCrossSectionIntegralPrecision_dX,
	int numberOfThreads
)
{
    //solve model at zero chemical potential up to some finite temperature
    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteTSolution = 
    solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(vacuum, minimumTemperature, numberOfPointsFromVacToMinTemp, gapPrecision, HYBRIDS);
    for (int i = 0; i < int(finiteTSolution.size()); ++i)
    {
        cout << finiteTSolution[i].getTemperature() << "\t" 
        	 << finiteTSolution[i].getUpQuarkEffectiveMass() << "\t"
        	 << finiteTSolution[i].getDownQuarkEffectiveMass() << "\t"
        	 << finiteTSolution[i].getStrangeQuarkEffectiveMass() << "\n";
    }


    vector<SU3NJL3DCutoffFixedChemPotTemp> finiteChemPotSolution = solveFromFiniteTemperatureToFiniteChemicalPotential(finiteTSolution[finiteTSolution.size()-1], chemPot, numberOfPointsMinTempToChemPot, gapPrecision, HYBRIDS);
    for (int i = 0; i < int(finiteChemPotSolution.size()); ++i)
    {   
        cout << finiteChemPotSolution[i].getTemperature() << "\t"
             << finiteChemPotSolution[i].getUpQuarkChemicalPotential() << "\t"
             << finiteChemPotSolution[i].getDownQuarkChemicalPotential() << "\t"
             << finiteChemPotSolution[i].getStrangeQuarkChemicalPotential() << "\t"
             << finiteChemPotSolution[i].getUpQuarkEffectiveMass() << "\t"
             << finiteChemPotSolution[i].getDownQuarkEffectiveMass() << "\t"
             << finiteChemPotSolution[i].getStrangeQuarkEffectiveMass() << "\n";
    }


    finiteTSolution = 
    solveFromLowToHighTemperature(finiteChemPotSolution[finiteChemPotSolution.size()-1], maximumTemperature, numberOfPointsFromMinToMaxTemp, gapPrecision, HYBRIDS);
    for (int i = 0; i < int(finiteTSolution.size()); ++i)
    {
        cout << finiteTSolution[i].getTemperature() << "\t" 
        	 << finiteTSolution[i].getUpQuarkChemicalPotential() << "\t"
        	 << finiteTSolution[i].getUpQuarkEffectiveMass() << "\t"
        	 << finiteTSolution[i].getDownQuarkEffectiveMass() << "\t"
        	 << finiteTSolution[i].getStrangeQuarkEffectiveMass() << "\n";
    }


    evaluateIsospinSymmetricIntegratedCrossSectionsAlongFixedChemicalPotentialTrajectory(
		finiteTSolution, 
        propagatorIntegralPrecision,
        largeAngleScatteringContribution, 
        crossSectionIntegralPrecision,
        integratedCrossSectionIntegralPrecision_dXdY, 
        integratedCrossSectionIntegralPrecision_dX, 
        approximationMethod,
		numberOfThreads
	);
}

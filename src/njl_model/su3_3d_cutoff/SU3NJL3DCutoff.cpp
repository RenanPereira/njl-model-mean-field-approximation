#include <cmath>
#include <iostream>

#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoff.h"
#include "physics_utils/distribution_functions.h"
#include "njl_model/n_fermion_line_integrals/one_fermion_line_integral_3d_cutoff.h"

using namespace std;

/*
Coupling constants at the Lagrangian Level
GSP4 = g/2;
Gdet = 8 kappa;

GSP81 = 16 g1;
GSP82 = 16 g2;

GVP4 = -gOmega;
GVIPI4 = -gRho;

GVP8 = -gOmega2;
GVIPI8 = -gRho2;
GVPVIPI8 = -gOmegaRho;

GSPVP8 = -gSigmaOmega;
GSPVIPI8 = -gSigmaRho;

GVP12 = -gOmega3;

GVP16 = -gOmega4;
*/

double SU3BaryonDensity(double upQuarkDensity, double downQuarkDensity, double strangeQuarkDensity)
{
	double rhoB = (1.0/3.0)*( upQuarkDensity + downQuarkDensity + strangeQuarkDensity );

	return rhoB;
}

double SU3NJLNulledGapEquation(
    NJLDimensionfulCouplings couplings, 
    double Mf_minus_m0f, 
    double sigma_f, 
    double sigma_fplus1, 
    double sigma_fplus2, 
    double rho_f, 
    double rho_fplus1, 
    double rho_fplus2
)
{	
    //no interactions
    double nullGap = 0.0;
    nullGap = Mf_minus_m0f;

    //4scalar interactions
    double gS = couplings.getFourQuarkSPCoupling();
    nullGap = nullGap 
            + 2.0*gS*sigma_f;

    //6scalar interactions
    double kappa = couplings.getDeterminantCoupling();
	nullGap = nullGap 
			+ 2.0*kappa*sigma_fplus1*sigma_fplus2;

    //8scalar interactions
    double g1 = couplings.getEightQuarkSPOziViolatingCoupling();
    double g2 = couplings.getEightQuarkSPNonOziViolatingCoupling();
    nullGap = nullGap 
            + 4.0*g1*sigma_f*( pow(sigma_f,2) + pow(sigma_fplus1,2) + pow(sigma_fplus2,2) ) 
            + 4.0*g2*pow(sigma_f,3);

    //4scalar_4vector interactions
    double gSigmaOmega = couplings.getEightQuarkSPVPCoupling();
    double gSigmaRho = couplings.getEightQuarkSPVIPICoupling();
    nullGap = nullGap 
            - (8.0/3.0)*gSigmaOmega*( pow(rho_f + rho_fplus1 + rho_fplus2, 2) )*sigma_f
            - (16.0/3.0)*gSigmaRho*( pow(rho_f,2) + pow(rho_fplus1,2) + pow(rho_fplus2,2) - rho_f*rho_fplus1 - rho_fplus1*rho_fplus2 - rho_f*rho_fplus2 )*sigma_f;

	return nullGap;
}

double SU3NJLQuarkChemicalPotential(
    NJLDimensionfulCouplings couplings, 
    double effectiveCP_f, 
    double rho_f, 
    double rho_fplus1, 
    double rho_fplus2, 
    double sigma_f, 
    double sigma_fplus1, 
    double sigma_fplus2
)
{	
    //no interactions
    double cp_f = 0.0;
    cp_f = effectiveCP_f;

    //4vector interactions
    double gOmega = couplings.getFourQuarkVPCoupling();
    double gRho = couplings.getFourQuarkVIPICoupling();
    cp_f = cp_f 
         + (4.0/3.0)*gOmega*(rho_f + rho_fplus1 + rho_fplus2)
         + (4.0/3.0)*gRho*(2.0*rho_f - rho_fplus1 - rho_fplus2);

    //8vector interactions
    double gOmega2 = couplings.getEightQuarkVPCoupling();
    double gRho2 = couplings.getEightQuarkVIPICoupling();
    double gOmegaRho = couplings.getEightQuarkVPVIPICoupling();
	cp_f = cp_f 
		 + (16.0/9.0)*gOmega2*pow( (rho_f + rho_fplus1 + rho_fplus2) , 3)
		 + (32.0/9.0)*gRho2*pow( (2.0*rho_f - rho_fplus1 - rho_fplus2) , 3)
		 + (32.0/3.0)*gRho2*(rho_f - rho_fplus1)*(rho_fplus2 - rho_f)*(2.0*rho_f - rho_fplus1 - rho_fplus2)
		 + (8.0/9.0)*gOmegaRho*(rho_f + rho_fplus1 + rho_fplus2)*(4.0*pow(rho_f,2) + pow(rho_fplus1,2) + pow(rho_fplus2,2) - rho_f*rho_fplus1 - rho_f*rho_fplus2 - 4.0*rho_fplus1*rho_fplus2);

    //4scalar_4vector interactions
    double gSigmaOmega = couplings.getEightQuarkSPVPCoupling();
    double gSigmaRho = couplings.getEightQuarkSPVIPICoupling();
    cp_f = cp_f 
         + (8.0/3.0)*gSigmaOmega*(rho_f + rho_fplus1 + rho_fplus2)*( pow(sigma_f,2) + pow(sigma_fplus1,2) + pow(sigma_fplus2,2) )
         + (8.0/3.0)*gSigmaRho*(2.0*rho_f - rho_fplus1 - rho_fplus2)*( pow(sigma_f,2) + pow(sigma_fplus1,2) + pow(sigma_fplus2,2) );

    //12vector interactions
    double gOmega3 = couplings.getTwelveQuarkVPCoupling();
    cp_f = cp_f 
         + (16.0/9.0)*gOmega3*pow( (rho_f + rho_fplus1 + rho_fplus2) , 5);

    //16vector interactions
    double gOmega4 = couplings.getSixteenQuarkVPCoupling();
    cp_f = cp_f 
         + (128.0/81.0)*gOmega4*pow( (rho_f + rho_fplus1 + rho_fplus2) , 7);

    //multi VP quark interaction     
    if ( couplings.getInteractionsIncludeMultiQuarkVPCouplings()==true )
    {
        double rho = rho_f + rho_fplus1 + rho_fplus2;
        for (int i = 0; i < couplings.numberOfMultiQuarkVPCoupling(); ++i)
        {   
            double X = i + 1.0;
            double gOmegaX = couplings.getMultiQuarkVPCoupling(i);
            cp_f = cp_f 
                 + ( 2.0*X*pow(2.0/3.0, X) )*gOmegaX*pow(rho, 2.0*X - 1.0);
        }
    }

	return cp_f;
}

double SU3NJLInteractionPotential(
    NJLDimensionfulCouplings couplings, 
    double sigmaU, 
    double sigmaD, 
    double sigmaS, 
    double rhoU, 
    double rhoD, 
    double rhoS
)
{	
    //no interactions
	double interactionPotential = 0.0;

    //4scalar interactions
    double gS = couplings.getFourQuarkSPCoupling();
    interactionPotential = interactionPotential 
                         + 1.0*gS*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) );

    //6scalar interaction
    double kappa = couplings.getDeterminantCoupling();
    interactionPotential = interactionPotential
                         + 4.0*kappa*sigmaU*sigmaD*sigmaS;

    //8scalar interactions
    double g1 = couplings.getEightQuarkSPOziViolatingCoupling();
    double g2 = couplings.getEightQuarkSPNonOziViolatingCoupling();
    interactionPotential = interactionPotential
                         + 3.0*g1*pow( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) ,2)
                         + 3.0*g2*( pow(sigmaU,4) + pow(sigmaD,4) + pow(sigmaS,4) );

    //4vector interactions
    double gOmega = couplings.getFourQuarkVPCoupling();
    double gRho = couplings.getFourQuarkVIPICoupling();
    interactionPotential = interactionPotential
                         - (2.0/3.0)*gOmega*pow(rhoU+rhoD+rhoS, 2)
                         - gRho*( pow(rhoU-rhoD, 2) + (1.0/3.0)*pow(rhoU+rhoD-2.0*rhoS, 2) );

    //8vector interactions
    double gOmega2 = couplings.getEightQuarkVPCoupling();
    double gRho2 = couplings.getEightQuarkVIPICoupling();
    double gOmegaRho = couplings.getEightQuarkVPVIPICoupling();
    interactionPotential = interactionPotential
                         - (4.0/3.0)*gOmega2*pow(rhoU+rhoD+rhoS, 4)
                         - 3.0*gRho2*pow( pow(rhoU-rhoD, 2) + (1.0/3.0)*pow(rhoU+rhoD-2.0*rhoS, 2) , 2)
                         - 2.0*gOmegaRho*pow(rhoU+rhoD+rhoS, 2)*( pow(rhoU-rhoD, 2) + (1.0/3.0)*pow(rhoU+rhoD-2.0*rhoS, 2) );

    //4scalar_4vector interactions
    double gSigmaOmega = couplings.getEightQuarkSPVPCoupling();
    double gSigmaRho = couplings.getEightQuarkSPVIPICoupling();
    interactionPotential = interactionPotential
                  		 - 4.0*gSigmaOmega*( pow(rhoU+rhoD+rhoS, 2) )*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
                  		 - 8.0*gSigmaRho*(  pow(rhoU,2) + pow(rhoD,2) + pow(rhoS,2) - rhoU*rhoD - rhoD*rhoS - rhoS*rhoU )*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) );

    //12vector interactions
    double gOmega3 = couplings.getTwelveQuarkVPCoupling();
    interactionPotential = interactionPotential 
                         - (40.0/27.0)*gOmega3*pow(rhoU+rhoD+rhoS, 6);       

    //16vector interactions
    double gOmega4 = couplings.getSixteenQuarkVPCoupling();
    interactionPotential = interactionPotential 
                         - (112.0/81.0)*gOmega4*pow(rhoU+rhoD+rhoS, 8);     

    //multi VP quark interaction     
    if ( couplings.getInteractionsIncludeMultiQuarkVPCouplings()==true )
    {
        double rho = rhoU+rhoD+rhoS;
        for (int i = 0; i < couplings.numberOfMultiQuarkVPCoupling(); ++i)
        {   
            double X = i + 1.0;
            double gOmegaX = couplings.getMultiQuarkVPCoupling(i);
            interactionPotential = interactionPotential 
                                 - ( (2.0*X - 1.0)*pow(2.0/3.0, X) )*gOmegaX*pow(rho, 2.0*X);
        }
    }

    return interactionPotential;
}

double fermionParticleDensityIntegrand(double k, void *parameters)
{   
    ThermodynamicsIntegrandParameters aux(parameters);
    double T = aux.getTemperature();
    double effCP = aux.getEffectiveChemicalPotential();
    double M = aux.getEffectiveMass();

    double E = sqrt( pow(k,2) + pow(M,2) );

    double integrand = pow(k,2)*( fermiDistribution(T, E-effCP) - fermiDistribution(T, E+effCP) );

    return integrand;
}

//CTmu term that correctly reproduces the Stefan–Boltzmann limit (particle density)
double fermionParticleDensity3DCutoffStefanBoltzmannCTmu(
    double cutoff, 
    double T, 
    double effChemPot, 
    double integralPrecision
)
{   
    int integralWorkspace = 1000;
    ThermodynamicsIntegrandParameters params("fermionParticleDensityCTmuIntegral", T, effChemPot, 0.0);

    double particleDensity = 0.0;
    if ( T>0.0 )
    {   
    	//finite T
        Integration1DimGSLQAGIU fermionParticleDensityCTmu(
            cutoff, 
            &params, 
            fermionParticleDensityIntegrand, 
            integralPrecision, integralPrecision, integralWorkspace
        );
        particleDensity = fermionParticleDensityCTmu.evaluate();
    }
    else
    {   
        //T=0 limit
        particleDensity = 0.0;
    }
    particleDensity = ( 1.0/(M_PI*M_PI) )*particleDensity;

    return particleDensity;
}

double fermionParticleDensity3DCutoff(
    NJL3DCutoffRegularizationScheme reguScheme, 
    double cutoff, 
    double T, 
    double effChemPot, 
    double effMass, 
    double integralPrecision
)
{   
    int integralWorkspace = 1000;
    ThermodynamicsIntegrandParameters params("fermionParticleDensity", T, effChemPot, effMass);

    double fermionDensity = 0.0;
    if ( T>0.0 )
    {   
        //finite T
        if ( reguScheme==CUTOFF_EVERYWHERE || reguScheme==CUTOFF_EVERYWHERE_WITH_CTMU )
        {   
            Integration1DimGSLQAGS fermionParticleDensity(
                0.0, 
                cutoff, 
                &params, 
                fermionParticleDensityIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionDensity = ( 1.0/(M_PI*M_PI) )*fermionParticleDensity.evaluate();
        }
        else if ( reguScheme==CUTOFF_ON_DIVERGENT_INTEGRALS_ONLY )
        {
            Integration1DimGSLQAGIU fermionParticleDensity(
                0.0, 
                &params, 
                fermionParticleDensityIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionDensity = ( 1.0/(M_PI*M_PI) )*fermionParticleDensity.evaluate();
        }
    }
    else
    {   
        //T=0 limit
        fermionDensity = pow(fermiMomentum(effChemPot,effMass),3)/( 3.0*M_PI*M_PI );
    }

    //if the chosen regularization includes the CTmu term, add it to the fermion density
	if ( reguScheme==CUTOFF_EVERYWHERE_WITH_CTMU )
	{
		fermionDensity = fermionDensity + fermionParticleDensity3DCutoffStefanBoltzmannCTmu(cutoff, T, effChemPot, integralPrecision);
	}

    return fermionDensity;
}

//Divergent contribution to the fermion pressure
double fermionPressureDivergentPrimitive3DCutoff(double k, double M)
{   
	double primitive = 0.0;
	
    double E = sqrt( pow(k,2) + pow(M,2) );
    
    primitive = (1.0/8.0)*( k*E*( 2.0*pow(k,2) + pow(M,2) ) - 0.5*pow(M,4)*( log( 1.0 + k/E ) - log( 1.0 - k/E ) ) );
	
    return primitive;
}

//Convergent contribution to the fermion pressure
double fermionPressureConvergentIntegrand(double k, void *parameters)
{   
    ThermodynamicsIntegrandParameters aux(parameters);
    double T = aux.getTemperature();
    double effCP = aux.getEffectiveChemicalPotential();
    double M = aux.getEffectiveMass();

    double E = sqrt( pow(k,2) + pow(M,2) );

    double integrand = pow(k,2)*( log( 1.0 + exp( -( E + effCP )/T ) ) + log( 1.0 + exp( -( E - effCP )/T ) ) );

    return integrand;
}

//CTmu term that correctly reproduces the Stefan–Boltzmann limit (pressure)
double fermionPressure3DCutoffStefanBoltzmannCTmu(double cutoff, double T, double effChemPot, double integralPrecision)
{   
    int integralWorkspace = 1000;
    ThermodynamicsIntegrandParameters params("fermionPressureCTmuIntegral", T, effChemPot, 0.0);

    double quarkPressure = 0.0;
    if ( T>0.0 )
    {   
    	//finite T
        Integration1DimGSLQAGIU fermionPressureCTmu(
            cutoff, 
            &params, 
            fermionPressureConvergentIntegrand, 
            integralPrecision, integralPrecision, integralWorkspace
        );
        quarkPressure = T*fermionPressureCTmu.evaluate();
    }
    else
    {   
        //T=0 limit
        quarkPressure = 0.0;
    }
    quarkPressure = ( 1.0/(M_PI*M_PI) )*quarkPressure;

    return quarkPressure;
}

double fermionPressure3DCutoff(
    NJL3DCutoffRegularizationScheme reguScheme, 
    double cutoff, 
    double T, 
    double effMass, 
    double effChemPot, 
    double integralPrecision
)
{   
    int integralWorkspace = 1000;
    ThermodynamicsIntegrandParameters params("fermionPressureConvergentIntegral", T, effChemPot, effMass);

    double fermionPressure = 0.0;
    if ( T>0.0 )
    {   
    	//finite T
    	if ( reguScheme==CUTOFF_EVERYWHERE || reguScheme==CUTOFF_EVERYWHERE_WITH_CTMU )
    	{  
            //divergent contribution
            fermionPressure = fermionPressureDivergentPrimitive3DCutoff(cutoff, effMass);

            //convergent contribution
            Integration1DimGSLQAGS fermionPressureConvergent(
                0.0, 
                cutoff, 
                &params, 
                fermionPressureConvergentIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionPressure = fermionPressure + T*fermionPressureConvergent.evaluate();
    	}
    	else if ( reguScheme==CUTOFF_ON_DIVERGENT_INTEGRALS_ONLY )
    	{  
            //divergent contribution
            fermionPressure = fermionPressureDivergentPrimitive3DCutoff(cutoff, effMass);

            //convergent contribution
            Integration1DimGSLQAGIU fermionPressureConvergent(
                0.0, 
                &params, 
                fermionPressureConvergentIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionPressure = fermionPressure + T*fermionPressureConvergent.evaluate();
    	}
    }
    else
    {
        //T=0 limit
        double fermiMom = fermiMomentum(effChemPot, effMass);
        fermionPressure = fermionPressureDivergentPrimitive3DCutoff(cutoff, effMass) - fermionPressureDivergentPrimitive3DCutoff(fermiMom, effMass);
        fermionPressure = fermionPressure + effChemPot*pow(fermiMom,3)/3.0;
    }
    fermionPressure = ( 1.0/(M_PI*M_PI) )*fermionPressure;

    //if the chosen regularization includes the CTmu term, add it to the pressure
	if ( reguScheme==CUTOFF_EVERYWHERE_WITH_CTMU )
	{
		fermionPressure = fermionPressure + fermionPressure3DCutoffStefanBoltzmannCTmu(cutoff, T, effChemPot, integralPrecision);
	}

    return fermionPressure;
}

//Convergent contribution to the fermion energy density
double fermionEnergyDensityConvergentIntegrand(double k, void *parameters)
{   
    ThermodynamicsIntegrandParameters aux(parameters);
    double T = aux.getTemperature();
    double effCP = aux.getEffectiveChemicalPotential();
    double cP = aux.getChemicalPotential();
    double M = aux.getEffectiveMass();

    double E = sqrt( pow(k,2) + pow(M,2) );

    double integrand = pow(k,2)*( ( effCP - cP - E )*fermiDistribution(T, E-effCP) + ( cP - effCP - E )*fermiDistribution(T, E+effCP) ); 

    return integrand;
}

double fermionEnergyDensity3DCutoff(
    NJL3DCutoffRegularizationScheme reguScheme, 
    double cutoff, 
    double T, 
    double effMass, 
    double effChemPot, 
    double chemPot, 
    double integralPrecision
)
{   
    int integralWorkspace = 1000;
    ThermodynamicsIntegrandParameters params("fermionEnergyConvergentIntegral", T, effChemPot, chemPot, effMass);

    double fermionEnergy = 0.0;
    if ( T>0.0 )
    {   
        //finite T
        if ( reguScheme==CUTOFF_EVERYWHERE || reguScheme==CUTOFF_EVERYWHERE_WITH_CTMU )
        {   
            //divergent contribution
            fermionEnergy = fermionPressureDivergentPrimitive3DCutoff(cutoff, effMass);

            //convergent contribution
            Integration1DimGSLQAGS fermionEnergyConvergent(
                0.0, 
                cutoff, 
                &params, 
                fermionEnergyDensityConvergentIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionEnergy = fermionEnergy + fermionEnergyConvergent.evaluate();
        }
        else if ( reguScheme==CUTOFF_ON_DIVERGENT_INTEGRALS_ONLY )
        {   
            //divergent contribution
            fermionEnergy = fermionPressureDivergentPrimitive3DCutoff(cutoff, effMass);

            //convergent contribution
            Integration1DimGSLQAGIU fermionEnergyConvergent(
                0.0, 
                &params, 
                fermionEnergyDensityConvergentIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionEnergy = fermionEnergy + fermionEnergyConvergent.evaluate();
        }
    }
    else
    {   
        //T=0 limit
        double fermiMom = fermiMomentum(effChemPot, effMass);
        fermionEnergy = fermionPressureDivergentPrimitive3DCutoff(cutoff, effMass) - fermionPressureDivergentPrimitive3DCutoff(fermiMom, effMass);
        fermionEnergy = fermionEnergy + ( effChemPot - chemPot )*pow(fermiMom,3)/3.0;
    }
    fermionEnergy = - ( 1.0/(M_PI*M_PI) )*fermionEnergy;

	//if the chosen regularization includes the CTmu term, add it to the energy density
	if ( reguScheme==CUTOFF_EVERYWHERE_WITH_CTMU )
	{
		fermionEnergy = fermionEnergy + fermionEnergyDensity3DCutoffStefanBoltzmannCTmu(cutoff, T, effChemPot, chemPot, integralPrecision);
	}
    
    return fermionEnergy;
}

double fermionEnergyDensity3DCutoffStefanBoltzmannCTmu(
    double cutoff, 
    double T, 
    double effChemPot, 
    double chemPot, 
    double integralPrecision
)
{   
    int integralWorkspace = 1000;
    ThermodynamicsIntegrandParameters params("fermionEnergyCTmuIntegral", T, effChemPot, chemPot, 0.0);

    double energy = 0.0;
    if ( T>0.0 )
    {   
    	//finite T
        Integration1DimGSLQAGIU fermionEnergyCTmu(
            cutoff, 
            &params, 
            fermionEnergyDensityConvergentIntegrand, 
            integralPrecision, integralPrecision, integralWorkspace
        );
        energy = fermionEnergyCTmu.evaluate();
    }
    else
    {   
        //T=0 limit
        energy = 0.0;
    }
    energy = -( 1.0/(M_PI*M_PI) )*energy;

    return energy;
}

//Convergent contribution to the fermion entropy
double fermionEntropyDensityConvergentIntegrand(double k, void *parameters)
{   
    ThermodynamicsIntegrandParameters aux(parameters);
    double T = aux.getTemperature();
    double effCP = aux.getEffectiveChemicalPotential();
    double M = aux.getEffectiveMass();

    double E = sqrt( pow(k,2) + pow(M,2) );

    double integrand = pow(k,2)*( ( (E-effCP)/T )*fermiDistribution(T, E-effCP) + ( (E+effCP)/T )*fermiDistribution(T, E+effCP) ); 

    return integrand;
}

double fermionEntropyDensity3DCutoff(
    NJL3DCutoffRegularizationScheme reguScheme, 
    double cutoff, 
    double T, 
    double effMass, 
    double effChemPot, 
    double integralPrecision
)
{   
    int integralWorkspace = 1000;
    ThermodynamicsIntegrandParameters params("fermionEntropyIntegral", T, effChemPot, effMass);

    double fermionEntropy = 0.0;
    if ( T>0.0 )
    {   
    	//finite T
        if ( reguScheme==CUTOFF_EVERYWHERE || reguScheme==CUTOFF_EVERYWHERE_WITH_CTMU )
        {   
            //convergent contributions
            Integration1DimGSLQAGS fermionEntropyConvergent1(
                0.0, 
                cutoff, 
                &params, 
                fermionPressureConvergentIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionEntropy = fermionEntropyConvergent1.evaluate();

            Integration1DimGSLQAGS fermionEntropyConvergent2(
                0.0, 
                cutoff, 
                &params, 
                fermionEntropyDensityConvergentIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionEntropy = fermionEntropy + fermionEntropyConvergent2.evaluate();
        }
        else if ( reguScheme==CUTOFF_ON_DIVERGENT_INTEGRALS_ONLY )
        {   
            //convergent contributions
            Integration1DimGSLQAGIU fermionEntropyConvergent1(
                0.0, 
                &params, 
                fermionPressureConvergentIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionEntropy = fermionEntropyConvergent1.evaluate();

            Integration1DimGSLQAGIU fermionEntropyConvergent2(
                0.0, 
                &params, 
                fermionEntropyDensityConvergentIntegrand, 
                integralPrecision, integralPrecision, integralWorkspace
            );
            fermionEntropy = fermionEntropy + fermionEntropyConvergent2.evaluate();
        }
    }
    else
    {   
    	//T=0 limit
        fermionEntropy = 0;
    }
    fermionEntropy = ( 1.0/(M_PI*M_PI) )*fermionEntropy;
    

	//if the chosen regularization includes the CTmu term, add it to the entropy density
	if ( reguScheme==CUTOFF_EVERYWHERE_WITH_CTMU )
	{
		fermionEntropy = fermionEntropy + fermionEntropyDensity3DCutoffStefanBoltzmannCTmu(cutoff, T, effChemPot, integralPrecision);
	}

    return fermionEntropy;
}

double fermionEntropyDensity3DCutoffStefanBoltzmannCTmu(double cutoff, double T, double effChemPot, double integralPrecision)
{   
	int integralWorkspace = 1000;
    ThermodynamicsIntegrandParameters params("fermionEntropyCTmuIntegral", T, effChemPot, 0.0);

    double entropy = 0.0;
    if ( T>0.0 )
    {   
    	//finite T
        Integration1DimGSLQAGIU fermionEntropyConvergent1(
            cutoff, 
            &params, 
            fermionPressureConvergentIntegrand, 
            integralPrecision, integralPrecision, integralWorkspace
        );
        entropy = fermionEntropyConvergent1.evaluate();

        Integration1DimGSLQAGIU fermionEntropyConvergent2(
            cutoff, 
            &params, 
            fermionEntropyDensityConvergentIntegrand, 
            integralPrecision, integralPrecision, integralWorkspace
        );
        entropy = entropy + fermionEntropyConvergent2.evaluate();
    }
    else
    {   
        //T=0 limit
        entropy = 0.0;
    }
    entropy = ( 1.0/(M_PI*M_PI) )*entropy;

    return entropy;
}

double SU3NJL3DCutoffPressure(
    SU3NJL3DCutoffParameters parametersNJL, 
    double T, 
    double effMassU, 
    double effMassD, 
    double effMassS, 
    double effChemPotU, 
    double effChemPotD, 
    double effChemPotS
)
{   
	NJLDimensionfulCouplings couplings = parametersNJL.getDimensionfulCouplings();

	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();

    double Nc = parametersNJL.getNumberOfColours();

    double sigmaIntegralPrecision = parametersNJL.getSigmaIntegralPrecision();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

    //Calculate the sigma field and density for each quark flavour
    double sigmaU = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effChemPotU, effMassU, sigmaIntegralPrecision);
    double sigmaD = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effChemPotD, effMassD, sigmaIntegralPrecision);
    double sigmaS = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effChemPotS, effMassS, sigmaIntegralPrecision);

	double rhoU = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effChemPotU, effMassU, thermoIntegralPrecision);
	double rhoD = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effChemPotD, effMassD, thermoIntegralPrecision);
	double rhoS = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effChemPotS, effMassS, thermoIntegralPrecision);

	//calculate integral contributions to the pressure for each quark flavour
	double fermionPressureU = Nc*fermionPressure3DCutoff(reguScheme, cutoff, T, effMassU, effChemPotU, thermoIntegralPrecision);
	double fermionPressureD = Nc*fermionPressure3DCutoff(reguScheme, cutoff, T, effMassD, effChemPotD, thermoIntegralPrecision);
	double fermionPressureS = Nc*fermionPressure3DCutoff(reguScheme, cutoff, T, effMassS, effChemPotS, thermoIntegralPrecision);

	//calculate the SU3 NJL model pressure
	double pressure = 0.0;
	pressure = - SU3NJLInteractionPotential(couplings, sigmaU, sigmaD, sigmaS, rhoU, rhoD, rhoS)
			   + fermionPressureU 
			   + fermionPressureD 
			   + fermionPressureS;

    return pressure;
}

double SU3NJL3DCutoffEnergyDensity(
    SU3NJL3DCutoffParameters parametersNJL, 
    double T, 
    double effMassU, 
    double effMassD, 
    double effMassS, 
    double effChemPotU, 
    double effChemPotD, 
    double effChemPotS
)
{   
    NJLDimensionfulCouplings couplings = parametersNJL.getDimensionfulCouplings();

    NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();

    double Nc = parametersNJL.getNumberOfColours();

    double sigmaIntegralPrecision = parametersNJL.getSigmaIntegralPrecision();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();
    
    //Calculate the sigma field and density for each quark flavour
    double sigmaU = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effChemPotU, effMassU, sigmaIntegralPrecision);
    double sigmaD = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effChemPotD, effMassD, sigmaIntegralPrecision);
    double sigmaS = sigmaNJL3DCutoff(reguScheme, cutoff, Nc, T, effChemPotS, effMassS, sigmaIntegralPrecision);

    double rhoU = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effChemPotU, effMassU, thermoIntegralPrecision);
    double rhoD = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effChemPotD, effMassD, thermoIntegralPrecision);
    double rhoS = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effChemPotS, effMassS, thermoIntegralPrecision);

	//calculate quark chemical potential from the effective ones
    double chemPotU = SU3NJLQuarkChemicalPotential(couplings, effChemPotU, rhoU, rhoD, rhoS, sigmaU, sigmaD, sigmaS);
    double chemPotD = SU3NJLQuarkChemicalPotential(couplings, effChemPotD, rhoD, rhoS, rhoU, sigmaD, sigmaS, sigmaU);
    double chemPotS = SU3NJLQuarkChemicalPotential(couplings, effChemPotS, rhoS, rhoU, rhoD, sigmaS, sigmaU, sigmaD);

	//calculate integral contributions to the energy density for each quark flavour
	double fermionEnergyU = Nc*fermionEnergyDensity3DCutoff(reguScheme, cutoff, T, effMassU, effChemPotU, chemPotU, thermoIntegralPrecision);
	double fermionEnergyD = Nc*fermionEnergyDensity3DCutoff(reguScheme, cutoff, T, effMassD, effChemPotD, chemPotD, thermoIntegralPrecision);
	double fermionEnergyS = Nc*fermionEnergyDensity3DCutoff(reguScheme, cutoff, T, effMassS, effChemPotS, chemPotS, thermoIntegralPrecision);

	//calculate the SU3 NJL model energy density
	double energy = 0.0;
	energy = + SU3NJLInteractionPotential(couplings, sigmaU, sigmaD, sigmaS, rhoU, rhoD, rhoS)
			 + fermionEnergyU 
			 + fermionEnergyD 
			 + fermionEnergyS;

    return energy;
}


double SU3NJL3DCutoffEntropyDensity(
    SU3NJL3DCutoffParameters parametersNJL, 
    double T, 
    double effMassU, 
    double effMassD, 
    double effMassS, 
    double effChemPotU, 
    double effChemPotD, 
    double effChemPotS
)
{   
	double Nc = parametersNJL.getNumberOfColours();

	NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

	//calculate integral contributions to the entropy density for each quark flavour
	double fermionEntropyU = Nc*fermionEntropyDensity3DCutoff(reguScheme, cutoff, T, effMassU, effChemPotU, thermoIntegralPrecision);
	double fermionEntropyD = Nc*fermionEntropyDensity3DCutoff(reguScheme, cutoff, T, effMassD, effChemPotD, thermoIntegralPrecision);
	double fermionEntropyS = Nc*fermionEntropyDensity3DCutoff(reguScheme, cutoff, T, effMassS, effChemPotS, thermoIntegralPrecision);

	//calculate the SU3 NJL model entropy density
	double entropy = 0.0;
	entropy = + fermionEntropyU 
			  + fermionEntropyD 
			  + fermionEntropyS;

    return entropy;
}

double SU3NJL3DCutoffQuarkFlavourDensity(
    SU3NJL3DCutoffParameters& parametersNJL, 
    double T, 
    double effMass, 
    double effChemPot
)
{
    NJL3DCutoffRegularizationScheme reguScheme = parametersNJL.getNJL3DCutoffRegularizationScheme();
    double cutoff = parametersNJL.getThreeMomentumCutoff();
    double Nc = parametersNJL.getNumberOfColours();    
    double thermoIntegralPrecision = parametersNJL.getThermoIntegralPrecision();

    double rho = Nc*fermionParticleDensity3DCutoff(reguScheme, cutoff, T, effChemPot, effMass, thermoIntegralPrecision);

    return rho;
}

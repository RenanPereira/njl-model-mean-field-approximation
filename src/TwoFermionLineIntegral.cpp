#include <cmath>
#include <gsl/gsl_complex_math.h>
#include "TwoFermionLineIntegral.h"


//Zero variables necessary in this file: TFLI=TwoFermionLineIntegral
const double TFLI_ZERO = 1E-12;
const double TFLI_ZERO_MASS_DIFFERENCE = 1E-8;


//Choose if the imaginary part of the f1 loop function is symmetric or antisymmetric (w.r.t. to omega-zero component of the external momentum)
//Symmetric - the shift is made in the masses as: M^2 -> M^2 - i*epsilon
//Antisymmetric - the shift is made in the pole (external momemtum) as: eta*omega -> eta*( omega + i*epsilon )
//If true the symmetric is used, if false, the antisymmetric is used
bool symmetricImaginaryPart = true;


//Function that selects the sign in front of the pair-creation/annihilation imaginary contributions depending on how the poles are shifted
double imag16Pi2f1PairSign(double eta, bool symmetric)
{
    double imagf1PairSign = 0.0;
    if (symmetric){ imagf1PairSign = +1.0; }
    else{ imagf1PairSign = -eta; }

    return imagf1PairSign;
}


//Function that selects the sign in front of the scattering imaginary contributions depending on how the poles are shifted
double imag16Pi2f1ScatSign(double eta, double epsilon, bool symmetric)
{   
    double imagf1ScatSign = 0.0;
    if (symmetric){ imagf1ScatSign = +1*( -sign(epsilon) ); }
    else{ imagf1ScatSign = -eta; }

    return imagf1ScatSign;
}



////////////////////////////////////////////////////////////////////////////////////////
//Functions for the integration region dEdepsilon


double dEdepsilon_epsilonMinPlus(double M1, double M2, double k, double E)
{   
    double epsilonMinPlus = 0.0;

    if ( fabs(2.0*E-k)>TFLI_ZERO && fabs(2.0*E+k)>TFLI_ZERO )
    {
        if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
        {
            double sqrtArg = pow(k,2)*( 4.0*pow(E,2) - pow(M1-M2,2) - pow(k,2) )*( 4.0*pow(E,2) - pow(M1+M2,2) - pow(k,2) );
    
            epsilonMinPlus = ( -2.0*E*( pow(M1,2) - pow(M2,2) ) + sqrt( fabs(sqrtArg) ) )/( 4.0*pow(E,2) - pow(k,2) );

            if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
            { 
                cout << "Argument of the square root in dEdepsilon_epsilonMinPlus is negative and larger then variable TFLI_ZERO!\n"; 
            }
        }
        else
        {   
            double sqrtArg = pow(k,2)*( 4.0*pow(E,2) - pow(k,2) )*( 4.0*pow(E,2) - 4.0*pow(M1,2) - pow(k,2) );

            epsilonMinPlus = sqrt( fabs(sqrtArg) )/( 4.0*pow(E,2) - pow(k,2) );

            if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
            { 
                cout << "Argument of the square root in dEdepsilon_epsilonMinPlus is negative and larger then variable TFLI_ZERO!\n"; 
            }
        }
    }
    else if( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
    {
        epsilonMinPlus = -( M1-M2 )*( M1+M2 )/( 4.0*E ) - ( 2.0*E*( pow(M1,2) + pow(M2,2) ) )/( ( M1-M2 )*( M1+M2 ) );
    }
    else
    {
        cout << "The function dEdepsilon_epsilonMinPlus is not defined for these inputs!\n";
    }

    return epsilonMinPlus;
}


double dEdepsilon_epsilonMinMinus(double M1, double M2, double k, double E)
{   
    double epsilonMinMinus = 0.0;

    if ( fabs(2.0*E-k)>TFLI_ZERO && fabs(2.0*E+k)>TFLI_ZERO )
    {
        if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
        {
            double sqrtArg = pow(k,2)*( 4.0*pow(E,2) - pow(M1-M2,2) - pow(k,2) )*( 4.0*pow(E,2) - pow(M1+M2,2) - pow(k,2) );

            epsilonMinMinus = ( -2.0*E*( pow(M1,2) - pow(M2,2) ) - sqrt( fabs(sqrtArg) ) )/( 4.0*pow(E,2) - pow(k,2) );

            if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
            { 
                cout << "Argument of the square root in dEdepsilon_epsilonMinMinus is negative and larger then variable TFLI_ZERO!\n"; 
            }
        }
        else
        {   
            double sqrtArg = pow(k,2)*( 4.0*pow(E,2) - pow(k,2) )*( 4.0*pow(E,2) - 4.0*pow(M1,2) - pow(k,2) );

            epsilonMinMinus = -sqrt( fabs(sqrtArg) )/( 4.0*pow(E,2) - pow(k,2) );

            if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
            { 
                cout << "Argument of the square root in dEdepsilon_epsilonMinMinus is negative and larger then variable TFLI_ZERO!\n"; 
            }
        }
    }
    else if( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
    {
        epsilonMinMinus = -( M1-M2 )*( M1+M2 )/( 4.0*E ) - ( 2.0*E*( pow(M1,2) + pow(M2,2) ) )/( ( M1-M2 )*( M1+M2 ) );
    }
    else
    {
        cout << "The function dEdepsilon_epsilonMinMinus is not defined for these inputs!\n";
    }

    return epsilonMinMinus;
}


double dEdepsilon_epsilonLambdaM1Plus(double cutoff, double M1, double E)
{
    double epsilonLambdaM1Plus = 2.0*E + 2.0*sqrt( pow(M1,2) + pow(cutoff,2) );

    return epsilonLambdaM1Plus;
}


double dEdepsilon_epsilonLambdaM1Minus(double cutoff, double M1, double E)
{
    double epsilonLambdaM1Minus = 2.0*E - 2.0*sqrt( pow(M1,2) + pow(cutoff,2) );

    return epsilonLambdaM1Minus;
}


double dEdepsilon_epsilonLambdaM2Plus(double cutoff, double M2, double E)
{
    double epsilonLambdaM2Plus = -2.0*E + 2.0*sqrt( pow(M2,2) + pow(cutoff,2) );

    return epsilonLambdaM2Plus;
}


double dEdepsilon_epsilonLambdaM2Minus(double cutoff, double M2, double E)
{
    double epsilonLambdaM2Minus = -2.0*E - 2.0*sqrt( pow(M2,2) + pow(cutoff,2) );

    return epsilonLambdaM2Minus;
}


double dEdepsilon_gamma1E(double cutoff, double M1, double M2, double k)
{
    double gamma1E = 0.5*( sqrt( pow(M1,2) + pow(k-cutoff,2) ) + sqrt( pow(M2,2) + pow(cutoff,2) ) );

    return gamma1E;
}


double dEdepsilon_gamma1epsilon(double cutoff, double M1, double M2, double k)
{
    double gamma1epsilon = -sqrt( pow(M1,2) + pow(k-cutoff,2) ) + sqrt( pow(M2,2) + pow(cutoff,2) );

    return gamma1epsilon;
}


double dEdepsilon_gamma2E(double cutoff, double M1, double M2, double k)
{
    double gamma2E = 0.5*( sqrt( pow(M2,2) + pow(k-cutoff,2) ) + sqrt( pow(M1,2) + pow(cutoff,2) ) );

    return gamma2E;
}


double dEdepsilon_gamma2epsilon(double cutoff, double M1, double M2, double k)
{
    double gamma2epsilon = sqrt( pow(M2,2) + pow(k-cutoff,2) ) - sqrt( pow(M1,2) + pow(cutoff,2) );

    return gamma2epsilon;
}


double dEdepsilon_alpha1E(double M1, double M2, double k)
{
    double alpha1E = 0.5*sqrt( pow(k,2) + pow(M1+M2,2) );

    return alpha1E;
}


double dEdepsilon_alpha1epsilon(double M1, double M2, double k)
{
    double alpha1epsilon = -( ( M1-M2 )*sqrt( pow(k,2) + pow(M1+M2,2) ) )/( M1+M2 );

    return alpha1epsilon;
}


double dEdepsilon_LambdaSwitchE(double cutoff, double M1, double M2)
{
    double LambdaSwitchE = 0.5*( sqrt( pow(M1,2) + pow(cutoff,2) ) + sqrt( pow(M2,2) + pow(cutoff,2) ) );

    return LambdaSwitchE;
}


double dEdepsilon_LambdaSwitchepsilon(double cutoff, double M1, double M2)
{
    double LambdaSwitchepsilon = -sqrt( pow(M1,2) + pow(cutoff,2) ) + sqrt( pow(M2,2) + pow(cutoff,2) );

    return LambdaSwitchepsilon;
}


double dEdepsilon_epsilonMax(double cutoff, double M1, double M2, double k, double E)
{
    double epsilonMax = 0.0;

    if ( E>dEdepsilon_gamma1E(cutoff, M1, M2, k) )
    {
        epsilonMax = dEdepsilon_epsilonLambdaM2Plus(cutoff, M2, E);
    }
    else
    {
        epsilonMax = dEdepsilon_epsilonMinPlus(M1, M2, k, E);
    }

    return epsilonMax;
}


double dEdepsilon_epsilonMin(double cutoff, double M1, double M2, double k, double E)
{
    double epsilonMin = 0.0;

    if ( E>dEdepsilon_gamma2E(cutoff, M1, M2, k) )
    {
        epsilonMin = dEdepsilon_epsilonLambdaM1Minus(cutoff, M1, E);
    }
    else
    {
        epsilonMin = dEdepsilon_epsilonMinMinus(M1, M2, k, E);
    }

    return epsilonMin;
}


double dEdepsilon_EMinM1LargerM2(double cutoff, double M1, double M2, double k)
{
    double EMin = 0.0;

    double alpha1E = dEdepsilon_alpha1E(M1, M2, k);
    double alpha1epsilon = dEdepsilon_alpha1epsilon(M1, M2, k);

    if ( dEdepsilon_epsilonLambdaM1Minus(cutoff, M1, alpha1E)>alpha1epsilon )
    {
        EMin = dEdepsilon_gamma2E(cutoff, M1, M2, k);
    }
    else
    {
        EMin = alpha1E;
    }

    return EMin;
}


double dEdepsilon_EMinM2LargerM1(double cutoff, double M1, double M2, double k)
{
    double EMin = 0.0;

    double alpha1E = dEdepsilon_alpha1E(M1, M2, k);
    double alpha1epsilon = dEdepsilon_alpha1epsilon(M1, M2, k);

    if ( dEdepsilon_epsilonLambdaM2Plus(cutoff, M2, alpha1E)<alpha1epsilon )
    {
        EMin = dEdepsilon_gamma1E(cutoff, M1, M2, k);
    }
    else
    {
        EMin = alpha1E;
    }

    return EMin;
}


double dEdepsilon_EMin(double cutoff, double M1, double M2, double k)
{
    double EMin = 0.0;

    if ( k>=2.0*cutoff )
    {
        EMin = 0.0;
    }
    else
    {
        if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
        {
            if (M1>M2)
            {
                EMin = dEdepsilon_EMinM1LargerM2(cutoff, M1, M2, k);
            }
            else
            {
                EMin = dEdepsilon_EMinM2LargerM1(cutoff, M1, M2, k);
            }
        }
        else
        {
            EMin = dEdepsilon_alpha1E(M1, M2, k);
        }
    }

    return EMin;
}


double dEdepsilon_EMax(double cutoff, double M1, double M2, double k)
{
    double EMax = 0.0;

    if ( k>=2.0*cutoff )
    {
        EMax = 0.0;
    }
    else
    {
        EMax = dEdepsilon_LambdaSwitchE(cutoff, M1, M2);
    }

    return EMax;
}


double dEdepsilon_AreaIntegrand(double E, void *parameters)
{
    TwoFermionLine3DCutoffIntegrand aux(parameters);
    double cutoff = aux.getThreeMomentumCutoff();
    double M1 = aux.getEffectiveMass1();
    double M2 = aux.getEffectiveMass2();
    double k = aux.getThreeMomentum();

    double AreaIntegrand = dEdepsilon_epsilonMax(cutoff, M1, M2, k, E) - dEdepsilon_epsilonMin(cutoff, M1, M2, k, E);

    return AreaIntegrand;
}


double dEdepsilon_Area(double cutoff, double M1, double M2, double k, double integralPrecision)
{   
    int integrationWorkspace = 1000;
    TwoFermionLine3DCutoffIntegrand params("dEdepsilon_Area", cutoff, M1, M2, k);

    double A = dEdepsilon_EMin(cutoff, M1, M2, k);
    double B = dEdepsilon_EMax(cutoff, M1, M2, k);

    Integration1DimGSLQAGS area(A, B, &params, dEdepsilon_AreaIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
    double dEdepsilonArea = area.evaluate();

    return dEdepsilonArea;
}


////////////////////////////////////////////////////////////////////////////////////////
//Functions for the integration region depsilondE


double depsilondE_EMinPlus(double M1, double M2, double k, double epsilon)
{   
    double EMinPlus = 0.0;

    if ( fabs(epsilon-k)>TFLI_ZERO && fabs(epsilon+k)>TFLI_ZERO )
    {
        if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
        {
            double sqrt_arg = pow(k,2)*( pow(k,2) + pow(M1-M2,2) - pow(epsilon,2) )*( pow(k,2) + pow(M1+M2,2) - pow(epsilon,2) );

            EMinPlus = ( epsilon*( pow(M1,2) - pow(M2,2) ) + sqrt( fabs(sqrt_arg) ) )/( 2.0*( k - epsilon )*( k + epsilon ) );

            if ( sqrt_arg<0 && fabs(sqrt_arg)>TFLI_ZERO )
            { 
                cout << "Argument of the square root in depsilondE_Eminplus is negative and larger then variable TFLI_ZERO!\n"; 
            }
        }
        else
        {   
            double sqrt_arg = pow(k,2)*( pow(k,2) - pow(epsilon,2) )*( pow(k,2) + 4.0*pow(M1,2) - pow(epsilon,2) );

            EMinPlus = 0.5*sqrt( sqrt_arg )/( pow(k,2) - pow(epsilon,2) );

            if ( sqrt_arg<0 && fabs(sqrt_arg)>TFLI_ZERO )
            { 
                cout << "Argument of the square root in depsilondE_Eminplus is negative and larger then variable TFLI_ZERO!\n"; 
            }
        }
    }
    else if( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
    {
        EMinPlus = - ( M1-M2 )*( M1+M2 )/( 4.0*epsilon ) - epsilon*( pow(M1,2) + pow(M2,2) )/( 2.0*( M1-M2 )*( M1+M2 ) );
    }
    else
    {
        cout << "The function depsilondE_Eminplus is not defined for these inputs!\n";
    }

    return EMinPlus;
}


double depsilondE_EMinMinus(double M1, double M2, double k, double epsilon)
{   
    double EMinMinus = 0.0;

    if ( fabs(epsilon-k)>TFLI_ZERO && fabs(epsilon+k)>TFLI_ZERO )
    {
        if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
        {
            double sqrt_arg = pow(k,2)*( pow(k,2) + pow(M1-M2,2) - pow(epsilon,2) )*( pow(k,2) + pow(M1+M2,2) - pow(epsilon,2) );

            EMinMinus = ( epsilon*( pow(M1,2) - pow(M2,2) ) - sqrt( fabs(sqrt_arg) ) )/( 2.0*( k-epsilon )*( k+epsilon ) );

            if ( sqrt_arg<0 && fabs(sqrt_arg)>TFLI_ZERO )
            { 
                cout << "Argument of the square root in depsilondE_Eminminus is negative and larger then variable TFLI_ZERO!\n"; 
            }
        }
        else
        {   
            double sqrt_arg = pow(k,2)*( pow(k,2) - pow(epsilon,2) )*( pow(k,2) + 4.0*pow(M1,2) - pow(epsilon,2) );

            EMinMinus = -0.5*sqrt( sqrt_arg )/( pow(k,2) - pow(epsilon,2) );

            if ( sqrt_arg<0 && fabs(sqrt_arg)>TFLI_ZERO )
            { 
                cout << "Argument of the square root in depsilondE_Eminminus is negative and larger then variable TFLI_ZERO!\n"; 
            }
        }
    }
    else if( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
    {
        EMinMinus = - ( M1-M2 )*( M1+M2 )/( 4.0*epsilon ) - epsilon*( pow(M1,2) + pow(M2,2) )/( 2.0*( M1-M2 )*( M1+M2 ) );
    }
    else
    {
        cout << "The function depsilondE_Eminminus is not defined for these inputs!\n";
    }

    return EMinMinus;
}


double depsilondE_ELambdaM1Plus(double cutoff, double M1, double epsilon)
{
    double ELambdaM1Plus = 0.5*epsilon + sqrt( pow(M1,2) + pow(cutoff,2) );

    return ELambdaM1Plus;
}


double depsilondE_ELambdaM1Minus(double cutoff, double M1, double epsilon)
{
    double ELambdaM1Minus = 0.5*epsilon - sqrt( pow(M1,2) + pow(cutoff,2) );

    return ELambdaM1Minus;
}


double depsilondE_ELambdaM2Plus(double cutoff, double M2, double epsilon)
{
    double ELambdaM2Plus = -0.5*epsilon + sqrt( pow(M2,2) + pow(cutoff,2) );

    return ELambdaM2Plus;
}


double depsilondE_ELambdaM2Minus(double cutoff, double M2, double epsilon)
{
    double ELambdaM2Minus = -0.5*epsilon - sqrt( pow(M2,2) + pow(cutoff,2) );

    return ELambdaM2Minus;
}


double depsilondE_LambdaSwitchepsilon(double cutoff, double M1, double M2)
{
    double LambdaSwitchepsilon = sqrt( pow(M2,2) + pow(cutoff,2) ) - sqrt( pow(M1,2) + pow(cutoff,2) );

    return LambdaSwitchepsilon;
}


double depsilondE_ELambdaPlus(double cutoff, double M1, double M2, double epsilon)
{   
    double ELambdaplus = 0.0;

    if ( epsilon>depsilondE_LambdaSwitchepsilon(cutoff, M1, M2) )
    {
        ELambdaplus = depsilondE_ELambdaM2Plus(cutoff, M2, epsilon);
    }
    else
    {
        ELambdaplus = depsilondE_ELambdaM1Plus(cutoff, M1, epsilon);
    }

    return ELambdaplus;
}


double depsilondE_gamma1epsilon(double cutoff, double M1, double M2, double k)
{
    double gamma1epsilon = sqrt( pow(M2,2) + pow(k-cutoff,2) ) - sqrt( pow(M1,2) + pow(cutoff,2) );

    return gamma1epsilon;
}


double depsilondE_gamma2epsilon(double cutoff, double M1, double M2, double k)
{
    double gamma2epsilon = -sqrt( pow(M1,2) + pow(k-cutoff,2) ) + sqrt( pow(M2,2) + pow(cutoff,2) );

    return gamma2epsilon;
}


double depsilondE_gamma1E(double cutoff, double M1, double M2, double k)
{
    double gamma1E = 0.5*( sqrt( pow(M2,2) + pow(k-cutoff,2) ) + sqrt( pow(M1,2) + pow(cutoff,2) ) );

    return gamma1E;
}


double depsilondE_gamma2E(double cutoff, double M1, double M2, double k)
{
    double gamma2E = 0.5*( sqrt( pow(M1,2) + pow(k-cutoff,2) ) + sqrt( pow(M2,2) + pow(cutoff,2) ) );

    return gamma2E;
}


double depsilondE_alpha1epsilon(double M1, double M2, double k)
{
    double alpha1epsilon = -sqrt( pow(M1-M2,2) + pow(k,2) );

    return alpha1epsilon;
}


double depsilondE_alpha2epsilon(double M1, double M2, double k)
{
    double alpha2epsilon = sqrt( pow(M1-M2,2) + pow(k,2) );

    return alpha2epsilon;
}


double depsilondE_alpha1E(double M1, double M2, double k)
{
    double alpha1E = 0.5*( ( M1+M2 )*sqrt( pow(M1-M2,2) + pow(k,2) ) )/( M1-M2 );

    return alpha1E;
}


double depsilondE_alpha2E(double M1, double M2, double k)
{
    double alpha2E = -0.5*( ( M1+M2 )*sqrt( pow(M1-M2,2) + pow(k,2) ) )/( M1-M2 );

    return alpha2E;
}


double depsilondE_EMaxM1LargerM2(double cutoff, double M1, double M2, double k, double epsilon)
{
    double EMax = 0.0;

    if ( epsilon<depsilondE_gamma1epsilon(cutoff, M1, M2, k) )
    {
        EMax = depsilondE_EMinMinus(M1, M2, k, epsilon);
    }
    else
    {
        EMax = depsilondE_ELambdaPlus(cutoff, M1, M2, epsilon);
    }

    return EMax;
}


double depsilondE_EMaxM2LargerM1(double cutoff, double M1, double M2, double k, double epsilon)
{
    double EMax = 0.0;

    if ( epsilon<depsilondE_gamma2epsilon(cutoff, M1, M2, k) )
    {
        EMax = depsilondE_ELambdaPlus(cutoff, M1, M2, epsilon);
    }
    else
    {
        EMax = depsilondE_EMinMinus(M1, M2, k, epsilon);
    }

    return EMax;
}


double depsilondE_EMax(double cutoff, double M1, double M2, double k, double epsilon)
{
    double EMax = 0.0;

    if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
    {
        if ( M1>M2 )
        {
            if ( depsilondE_alpha1E(M1, M2, k)<depsilondE_gamma1E(cutoff, M1, M2, k) )
            {
                EMax = depsilondE_EMaxM1LargerM2(cutoff, M1, M2, k, epsilon);
            }
            else
            {
                EMax = depsilondE_ELambdaPlus(cutoff, M1, M2, epsilon);
            }
        }
        else
        {
            if ( depsilondE_alpha2E(M1, M2, k)<depsilondE_gamma2E(cutoff, M1, M2, k) )
            {
                EMax = depsilondE_EMaxM2LargerM1(cutoff, M1, M2, k, epsilon);
            }
            else
            {
                EMax = depsilondE_ELambdaPlus(cutoff, M1, M2, epsilon);
            } 
        }
    }
    else
    {
        EMax = depsilondE_ELambdaPlus(cutoff, M1, M2, epsilon);
    }

    return EMax;
}


double depsilondE_EMin(double M1, double M2, double k, double epsilon)
{
    double EMin = 0.0;

    EMin = depsilondE_EMinPlus(M1, M2, k, epsilon);

    return EMin;
}


double depsilondE_epsilonMin(double cutoff, double M1, double M2, double k)
{
    double epsilonMin = 0.0;

    if ( k>=2.0*cutoff )
    {
        epsilonMin = 0.0;
    }
    else
    {
        if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
        {
            if (M1>M2)
            {   
                double alpha1epsilon = depsilondE_alpha1epsilon(M1, M2, k);

                if ( depsilondE_ELambdaPlus(cutoff, M1, M2, alpha1epsilon)>=depsilondE_alpha1E(M1, M2, k) )
                {
                    epsilonMin = alpha1epsilon;
                }
                else
                {
                    epsilonMin = depsilondE_gamma1epsilon(cutoff, M1, M2, k);
                }
            }
            else
            {
                epsilonMin = depsilondE_gamma1epsilon(cutoff, M1, M2, k);
            }
        }
        else
        {
            epsilonMin = depsilondE_gamma1epsilon(cutoff, M1, M2, k);
        }
    }

    return epsilonMin;
}


double depsilondE_epsilonMax(double cutoff, double M1, double M2, double k)
{
    double epsilonMax = 0.0;

    if ( k>=2.0*cutoff )
    {
        epsilonMax = 0.0;
    }
    else
    {
        if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
        {
            if (M2>M1)
            {   
                double alpha2epsilon = depsilondE_alpha2epsilon(M1, M2, k);

                if ( depsilondE_ELambdaPlus(cutoff, M1, M2, alpha2epsilon)>=depsilondE_alpha2E(M1, M2, k) )
                {
                    epsilonMax = alpha2epsilon;
                }
                else
                {
                    epsilonMax = depsilondE_gamma2epsilon(cutoff, M1, M2, k);
                }
            }
            else
            {
                epsilonMax = depsilondE_gamma2epsilon(cutoff, M1, M2, k);
            }
        }
        else
        {
            epsilonMax = depsilondE_gamma2epsilon(cutoff, M1, M2, k);
        }
    }

    return epsilonMax;
}


double depsilondE_AreaIntegrand(double epsilon, void *parameters)
{
    TwoFermionLine3DCutoffIntegrand aux(parameters);
    double cutoff = aux.getThreeMomentumCutoff();
    double M1 = aux.getEffectiveMass1();
    double M2 = aux.getEffectiveMass2();
    double k = aux.getThreeMomentum();

    double AreaIntegrand = depsilondE_EMax(cutoff, M1, M2, k, epsilon) - depsilondE_EMin(M1, M2, k, epsilon);

    return AreaIntegrand;
}


double depsilondE_Area(double cutoff, double M1, double M2, double k, double integralPrecision)
{
    int integrationWorkspace = 1000;
    TwoFermionLine3DCutoffIntegrand params("depsilondE_AREA", cutoff, M1, M2, k);

    double A = depsilondE_epsilonMin(cutoff, M1, M2, k);
    double B = depsilondE_epsilonMax(cutoff, M1, M2, k);

    Integration1DimGSLQAGS area(A, B, &params, depsilondE_AreaIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
    double depsilondEArea = area.evaluate();

    return depsilondEArea;
}


////////////////////////////////////////////////////////////////////////////////////////
// calculate real16Pi2f1Pair3DCutoff and imag16Pi2f1Pair3DCutoff using dEdepsilon


double gPlusEta(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double k, double eta, double E)
{   
    double C = dEdepsilon_epsilonMin(cutoff, M1, M2, k, E);
    double D = dEdepsilon_epsilonMax(cutoff, M1, M2, k, E);

/*
    double gPlusEtaAux;
    gPlusEtaAux =  + ( D - C )
                   - 2.0*T*log( 1 + exp( (- E + D/2 - eta*effCP1 )/T ) )
                   + 2.0*T*log( 1 + exp( (- E + C/2 - eta*effCP1 )/T ) )
                   + 2.0*T*log( 1 + exp( (- E - D/2 + eta*effCP2 )/T ) )
                   - 2.0*T*log( 1 + exp( (- E - C/2 + eta*effCP2 )/T ) );
*/

    //arguments inside the logarithms
    double argD1 = - E + 0.5*D - eta*effCP1;
    double argC1 = - E + 0.5*C - eta*effCP1;
    double argD2 = - E - 0.5*D + eta*effCP2;
    double argC2 = - E - 0.5*C + eta*effCP2;

    // + ( D - C )
    double gPlusEtaAux = ( D - C );
    
    // - 2.0*T*log( 1 + exp( (- E + D/2 - eta*effCP1 )/T ) )
    if ( isinf( exp(argD1/T) )==0 )
    { 
        gPlusEtaAux = gPlusEtaAux - 2.0*T*log( 1.0 + exp(argD1/T) ); 
    }
    else if( argD1>0 )
    { 
        gPlusEtaAux = gPlusEtaAux - 2.0*puiseuxExpansionTln1plusExpArgOverT(T, argD1); 
    }

    // + 2.0*T*log( 1 + exp( (- E + C/2 - eta*effCP1 )/T ) )
    if ( isinf( exp(argC1/T) )==0 )
    { 
        gPlusEtaAux = gPlusEtaAux + 2.0*T*log( 1.0 + exp(argC1/T) ); 
    }
    else if( argC1>0 )
    { 
        gPlusEtaAux = gPlusEtaAux + 2.0*puiseuxExpansionTln1plusExpArgOverT(T, argC1); 
    }

    // + 2.0*T*log( 1 + exp( (- E - D/2 + eta*effCP2 )/T ) )
    if ( isinf( exp(argD2/T) )==0 )
    { 
        gPlusEtaAux = gPlusEtaAux + 2.0*T*log( 1.0 + exp( argD2/T ) ); 
    }
    else if( argD2>0 )
    { 
        gPlusEtaAux = gPlusEtaAux + 2.0*puiseuxExpansionTln1plusExpArgOverT(T, argD2); 
    }

    // - 2.0*T*log( 1 + exp( (- E - C/2 + eta*effCP2 )/T ) )
    if ( isinf( exp(argC2/T) )==0 )
    { 
        gPlusEtaAux = gPlusEtaAux - 2.0*T*log( 1.0 + exp( argC2/T ) ); 
    }
    else if( argC2>0 )
    { 
        gPlusEtaAux = gPlusEtaAux - 2.0*puiseuxExpansionTln1plusExpArgOverT(T, argC2); 
    }

    return gPlusEtaAux;
}


double real16Pi2f1Pair3DCutoffNumerator(double E, void *parameters)
{   
    TwoFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double effCP1 = aux.getEffectiveChemicalPotential1();
    double effCP2 = aux.getEffectiveChemicalPotential2();
    double cutoff = aux.getThreeMomentumCutoff();
    double M1 = aux.getEffectiveMass1();
    double M2 = aux.getEffectiveMass2();
    double k = aux.getThreeMomentum();
    double eta = aux.getEtaVariable();

    double real16Pi2f1PairNumerator = gPlusEta(T, effCP1, effCP2, cutoff, M1, M2, k, eta, E);

    return real16Pi2f1PairNumerator;
}


double real16Pi2f1Pair3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double k, double integralPrecision)
{   
    int integrationWorkspace = 1000;

    double A = dEdepsilon_EMin(cutoff, M1, M2, k);
    double B = dEdepsilon_EMax(cutoff, M1, M2, k);
    double sing = 0.5*( w + (-effCP2 + effCP1) );

    double g1 = dEdepsilon_gamma1E(cutoff, M1, M2, k);
    double g2 = dEdepsilon_gamma2E(cutoff, M1, M2, k);

    TwoFermionLine3DCutoffIntegrand aux(T, effCP1, effCP2, cutoff, M1, M2, w, k);
    double real16Pi2f1Pair = 0.0;

    if ( M2>M1 && fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE ) //in this case g2>=g1
    {   
        ////////////////////////////////////////
        //eta = +1.0
        aux.setEtaVariable(+1.0);

        aux.setIntegralID("realf1PairPlus_A_g1");
        Integration1DimGSLQAWCQAGS realf1PairPlus_A_g1(A, g1, -sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairPlus_A_g1.evaluate();

        aux.setIntegralID("realf1PairPlus_g1_g2");
        Integration1DimGSLQAWCQAGS realf1PairPlus_g1_g2(g1, g2, -sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairPlus_g1_g2.evaluate();

        aux.setIntegralID("realf1PairPlus_g2_B");
        Integration1DimGSLQAWCQAGS realf1PairPlus_g2_B(g2, B, -sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairPlus_g2_B.evaluate();

        ////////////////////////////////////////
        //eta = -1.0
        aux.setEtaVariable(-1.0);

        aux.setIntegralID("realf1PairMinus_A_g1");
        Integration1DimGSLQAWCQAGS realf1PairMinus_A_g1(A, g1, +sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairMinus_A_g1.evaluate();

        aux.setIntegralID("realf1PairMinus_g1_g2");
        Integration1DimGSLQAWCQAGS realf1PairMinus_g1_g2(g1, g2, +sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairMinus_g1_g2.evaluate();

        aux.setIntegralID("realf1PairMinus_g2_B");
        Integration1DimGSLQAWCQAGS realf1PairMinus_g2_B(g2, B, +sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairMinus_g2_B.evaluate();
    }
    else if( M1>M2 && fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )  //in this case g1>=g2
    {
        ////////////////////////////////////////
        //eta = +1.0
        aux.setEtaVariable(+1.0);

        aux.setIntegralID("realf1PairPlus_A_g2");
        Integration1DimGSLQAWCQAGS realf1PairPlus_A_g2(A, g2, -sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairPlus_A_g2.evaluate();

        aux.setIntegralID("realf1PairPlus_g2_g1");
        Integration1DimGSLQAWCQAGS realf1PairPlus_g2_g1(g2, g1, -sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairPlus_g2_g1.evaluate();

        aux.setIntegralID("realf1PairPlus_g1_B");
        Integration1DimGSLQAWCQAGS realf1PairPlus_g1_B(g1, B, -sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairPlus_g1_B.evaluate();

        ////////////////////////////////////////
        //eta = -1.0
        aux.setEtaVariable(-1.0);

        aux.setIntegralID("realf1PairMinus_A_g2");
        Integration1DimGSLQAWCQAGS realf1PairMinus_A_g2(A, g2, +sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairMinus_A_g2.evaluate();

        aux.setIntegralID("realf1PairMinus_g2_g1");
        Integration1DimGSLQAWCQAGS realf1PairMinus_g2_g1(g2, g1, +sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairMinus_g2_g1.evaluate();

        aux.setIntegralID("realf1PairMinus_g1_B");
        Integration1DimGSLQAWCQAGS realf1PairMinus_g1_B(g1, B, +sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairMinus_g1_B.evaluate();
    }
    else //M1==M2, in this case g2==g1
    {
        ////////////////////////////////////////
        //eta = +1.0
        aux.setEtaVariable(+1.0);

        aux.setIntegralID("realf1PairPlus_A_g1");
        Integration1DimGSLQAWCQAGS realf1PairPlus_A_g1(A, g1, -sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairPlus_A_g1.evaluate();

        aux.setIntegralID("realf1PairPlus_g1_B");
        Integration1DimGSLQAWCQAGS realf1PairPlus_g1_B(g1, B, -sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairPlus_g1_B.evaluate();

        ////////////////////////////////////////
        //eta = -1.0
        aux.setEtaVariable(-1.0);

        aux.setIntegralID("realf1PairMinus_A_g1");
        Integration1DimGSLQAWCQAGS realf1PairMinus_A_g1(A, g1, +sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairMinus_A_g1.evaluate();

        aux.setIntegralID("realf1PairMinus_g1_B");
        Integration1DimGSLQAWCQAGS realf1PairMinus_g1_B(g1, B, +sing, &aux, real16Pi2f1Pair3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Pair = real16Pi2f1Pair + realf1PairMinus_g1_B.evaluate();

    }

    real16Pi2f1Pair = ( 1.0/(2.0*k) )*real16Pi2f1Pair;

    return real16Pi2f1Pair;
}


double imag16Pi2f1Pair3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double k)
{   
    double A = dEdepsilon_EMin(cutoff, M1, M2, k);
    double B = dEdepsilon_EMax(cutoff, M1, M2, k);
    double sing = 0.5*( w + (-effCP2 + effCP1) );

    double imag16Pi2f1Pair = 0.0;

    ////////////////////////////////////////
    //eta = +1.0
    double boxCarEtaPlus = heavisideTheta( -sing - A ) - heavisideTheta( -sing - B );
    if ( fabs(boxCarEtaPlus)>0 )
    {
        imag16Pi2f1Pair = imag16Pi2f1Pair + imag16Pi2f1PairSign(+1.0, symmetricImaginaryPart)*gPlusEta(T, effCP1, effCP2, cutoff, M1, M2, k, +1.0, -sing)*boxCarEtaPlus;
    }

    ////////////////////////////////////////
    //eta = -1.0
    double boxCarEtaMinus = heavisideTheta( sing - A ) - heavisideTheta( sing - B );
    if ( fabs(boxCarEtaMinus)>0 )
    {
        imag16Pi2f1Pair = imag16Pi2f1Pair + imag16Pi2f1PairSign(-1.0, symmetricImaginaryPart)*gPlusEta(T, effCP1, effCP2, cutoff, M1, M2, k, -1.0, +sing)*boxCarEtaMinus;
    }

    imag16Pi2f1Pair = ( M_PI/(2.0*k) )*imag16Pi2f1Pair;

    return imag16Pi2f1Pair;
}


////////////////////////////////////////////////////////////////////////////////////////
// real16Pi2f1Scat3DCutoff and imag16Pi2f1Scat3DCutoff using depsilondE


double gMinusEta(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double k, double eta, double epsilon)
{   
    double c = depsilondE_EMin(M1, M2, k, epsilon);
    double d = depsilondE_EMax(cutoff, M1, M2, k, epsilon);

/*
    double gMinusEtaAux;
    gMinusEtaAux = + T*log( 1 + exp( ( c - epsilon/2 - eta*effCP1 )/T ) )
                   - T*log( 1 + exp( ( d - epsilon/2 - eta*effCP1 )/T ) )
                   - T*log( 1 + exp( ( c + epsilon/2 - eta*effCP2 )/T ) )
                   + T*log( 1 + exp( ( d + epsilon/2 - eta*effCP2 )/T ) );
*/

    //arguments inside the logarithms
    double argD2 = d + 0.5*epsilon - eta*effCP2;
    double argD1 = d - 0.5*epsilon - eta*effCP1;
    double argC1 = c - 0.5*epsilon - eta*effCP1;
    double argC2 = c + 0.5*epsilon - eta*effCP2;

    double gMinusEtaAux = 0.0;

    // + T*log( 1 + exp( ( d + epsilon/2 - eta*effCP2 )/T ) )
    if ( isinf( exp(argD2/T) )==0 )
    { 
        gMinusEtaAux = gMinusEtaAux + T*log( 1 + exp(argD2/T) ); 
    }
    else if( argD2>0 )
    { 
        gMinusEtaAux = gMinusEtaAux + puiseuxExpansionTln1plusExpArgOverT(T, argD2); 
    }

    // - T*log( 1 + exp( ( d - epsilon/2 - eta*effCP1 )/T ) )
    if ( isinf( exp(argD1/T) )==0 )
    { 
        gMinusEtaAux = gMinusEtaAux - T*log( 1 + exp(argD1/T) ); 
    }
    else if( argD1>0 )
    { 
        gMinusEtaAux = gMinusEtaAux - puiseuxExpansionTln1plusExpArgOverT(T, argD1);
    }

    // + T*log( 1 + exp( ( c - epsilon/2 - eta*effCP1 )/T ) )
    if ( isinf( exp(argC1/T) )==0 )
    { 
        gMinusEtaAux = gMinusEtaAux + T*log( 1 + exp( argC1/T ) ); 
    }
    else if( argC1>0 )
    { 
        gMinusEtaAux = gMinusEtaAux + puiseuxExpansionTln1plusExpArgOverT(T, argC1);
    }

    // - T*log( 1 + exp( ( c + epsilon/2 - eta*effCP2 )/T ) )
    if ( isinf( exp(argC2/T) )==0 )
    { 
        gMinusEtaAux = gMinusEtaAux - T*log( 1 + exp( argC2/T ) ); 
    }
    else if( argC2>0 )
    { 
        gMinusEtaAux = gMinusEtaAux - puiseuxExpansionTln1plusExpArgOverT(T, argC2);
    }

    return gMinusEtaAux;
}


double real16Pi2f1Scat3DCutoffNumerator(double epsilon, void *parameters)
{   
    TwoFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double effCP1 = aux.getEffectiveChemicalPotential1();
    double effCP2 = aux.getEffectiveChemicalPotential2();
    double cutoff = aux.getThreeMomentumCutoff();
    double M1 = aux.getEffectiveMass1();
    double M2 = aux.getEffectiveMass2();
    double k = aux.getThreeMomentum();
    double eta = aux.getEtaVariable();

    double real16Pi2f1ScatNumerator = gMinusEta(T, effCP1, effCP2, cutoff, M1, M2, k, eta, epsilon);

    return real16Pi2f1ScatNumerator;
}


double real16Pi2f1Scat3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double k, double integralPrecision)
{   
    int integrationWorkspace = 1000;

    double a = depsilondE_epsilonMin(cutoff, M1, M2, k);
    double b = depsilondE_epsilonMax(cutoff, M1, M2, k);
    double sing = ( w + (-effCP2 + effCP1) );
    
    double g1 = depsilondE_gamma1epsilon(cutoff, M1, M2, k);
    double g2 = depsilondE_gamma2epsilon(cutoff, M1, M2, k);
    double L = depsilondE_LambdaSwitchepsilon(cutoff, M1, M2);

    TwoFermionLine3DCutoffIntegrand aux(T, effCP1, effCP2, cutoff, M1, M2, w, k);
    double real16Pi2f1Scat = 0.0;

    if ( M2>M1 && fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE ) //in this case g2>=g1
    {   
        ////////////////////////////////////////
        //eta = +1.0
        aux.setEtaVariable(+1.0);

        aux.setIntegralID("realf1ScatPlus_a_L");
        Integration1DimGSLQAWCQAGS realf1ScatPlus_a_L(a, L, -sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatPlus_a_L.evaluate();

        aux.setIntegralID("realf1ScatPlus_L_g2");
        Integration1DimGSLQAWCQAGS realf1ScatPlus_L_g2(L, g2, -sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatPlus_L_g2.evaluate();

        aux.setIntegralID("realf1ScatPlus_g2_b");
        Integration1DimGSLQAWCQAGS realf1ScatPlus_g2_b(g2, b, -sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatPlus_g2_b.evaluate();

        ////////////////////////////////////////
        //eta = -1.0
        aux.setEtaVariable(-1.0);

        aux.setIntegralID("realf1ScatMinus_a_L");
        Integration1DimGSLQAWCQAGS realf1ScatMinus_a_L(a, L, +sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatMinus_a_L.evaluate();

        aux.setIntegralID("realf1ScatMinus_L_g2");
        Integration1DimGSLQAWCQAGS realf1ScatMinus_L_g2(L, g2, +sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatMinus_L_g2.evaluate();

        aux.setIntegralID("realf1ScatMinus_g2_b");
        Integration1DimGSLQAWCQAGS realf1ScatMinus_g2_b(g2, b, +sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatMinus_g2_b.evaluate();
    }
    else if ( M1>M2 && fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )  //in this case, g1 is inside the interval of integration and g2 is not; L also separates the integration but g1<L;
    {
        ////////////////////////////////////////
        //eta = +1.0
        aux.setEtaVariable(+1.0);

        aux.setIntegralID("realf1ScatPlus_a_g1");
        Integration1DimGSLQAWCQAGS realf1ScatPlus_a_g1(a, g1, -sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatPlus_a_g1.evaluate();

        aux.setIntegralID("realf1ScatPlus_g1_L");
        Integration1DimGSLQAWCQAGS realf1ScatPlus_g1_L(g1, L, -sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatPlus_g1_L.evaluate();

        aux.setIntegralID("realf1ScatPlus_L_b");
        Integration1DimGSLQAWCQAGS realf1ScatPlus_L_b(L, b, -sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatPlus_L_b.evaluate();

        ////////////////////////////////////////
        //eta = -1.0
        aux.setEtaVariable(-1.0);

        aux.setIntegralID("realf1ScatMinus_a_g1");
        Integration1DimGSLQAWCQAGS realf1ScatMinus_a_g1(a, g1, +sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatMinus_a_g1.evaluate();

        aux.setIntegralID("realf1ScatMinus_g1_L");
        Integration1DimGSLQAWCQAGS realf1ScatMinus_g1_L(g1, L, +sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatMinus_g1_L.evaluate();

        aux.setIntegralID("realf1ScatMinus_L_b");
        Integration1DimGSLQAWCQAGS realf1ScatMinus_L_b(L, b, +sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatMinus_L_b.evaluate();
    }
    else  //M1==M2, in this case g1==a and g2=b, only the L point separates the integration
    {   
        ////////////////////////////////////////
        //eta = +1.0
        aux.setEtaVariable(+1.0);

        aux.setIntegralID("realf1ScatPlus_a_L");
        Integration1DimGSLQAWCQAGS realf1ScatPlus_a_L(a, L, -sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatPlus_a_L.evaluate();

        aux.setIntegralID("realf1ScatPlus_L_b");
        Integration1DimGSLQAWCQAGS realf1ScatPlus_L_b(L, b, -sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatPlus_L_b.evaluate();

        ////////////////////////////////////////
        //eta = -1.0
        aux.setEtaVariable(-1.0);

        aux.setIntegralID("realf1ScatMinus_a_L");
        Integration1DimGSLQAWCQAGS realf1ScatMinus_a_L(a, L, +sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatMinus_a_L.evaluate();

        aux.setIntegralID("realf1ScatMinus_L_b");
        Integration1DimGSLQAWCQAGS realf1ScatMinus_L_b(L, b, +sing, &aux, real16Pi2f1Scat3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace, TFLI_ZERO);
        real16Pi2f1Scat = real16Pi2f1Scat + realf1ScatMinus_L_b.evaluate();
    }

    real16Pi2f1Scat = -( 1.0/k )*real16Pi2f1Scat;

    return real16Pi2f1Scat;
}


double imag16Pi2f1Scat3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double k)
{
    double a = depsilondE_epsilonMin(cutoff, M1, M2, k);
    double b = depsilondE_epsilonMax(cutoff, M1, M2, k);
    double sing = ( w + (-effCP2 + effCP1) );

    double imag16Pi2f1Scat = 0.0;

    ////////////////////////////////////////
    //eta = +1.0
    double boxCarEtaPlus = heavisideTheta( -sing - a ) - heavisideTheta( -sing - b );
    if ( fabs(boxCarEtaPlus)>0 )
    {
        imag16Pi2f1Scat = imag16Pi2f1Scat + imag16Pi2f1ScatSign(+1.0, -sing, symmetricImaginaryPart)*gMinusEta(T, effCP1, effCP2, cutoff, M1, M2, k, +1.0, -sing)*boxCarEtaPlus;
    }

    ////////////////////////////////////////
    //eta = -1.0
    double boxCarEtaMinus = heavisideTheta( sing - a ) - heavisideTheta( sing - b );
    if ( fabs(boxCarEtaMinus)>0 )
    {
        imag16Pi2f1Scat = imag16Pi2f1Scat + imag16Pi2f1ScatSign(-1.0, +sing, symmetricImaginaryPart)*gMinusEta(T, effCP1, effCP2, cutoff, M1, M2, k, -1.0, +sing)*boxCarEtaMinus;
    }

    imag16Pi2f1Scat = -( M_PI/k )*imag16Pi2f1Scat;

    return imag16Pi2f1Scat;
}


////////////////////////////////////////////////////////////////////////////////////////
// real16Pi2f13DCutoffThreeMomentum and imag16Pi2f13DCutoffThreeMomentum


double real16Pi2f1Finite3Momentum3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double k, double integralPrecision)
{   
    double real16Pi2f13DCutoff = real16Pi2f1Pair3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, k, integralPrecision)
                               + real16Pi2f1Scat3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, k, integralPrecision);

    return real16Pi2f13DCutoff;
}


double imag16Pi2f1Finite3Momentum3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double k)
{   
    double imag16Pi2f13DCutoff = imag16Pi2f1Pair3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, k)
                               + imag16Pi2f1Scat3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, k);

    return imag16Pi2f13DCutoff;
}


////////////////////////////////////////////////////////////////////////////////////////
// calculate the special case in which the external three momentum is zero of the f1pair loop function
// lim k->0 f1pair: real16Pi2f1Pair3DCutoff and imag16Pi2f1Pair3DCutoff


double E0(double M1, double M2)
{
    double E0Aux = 0.5*( M1+M2 );

    return E0Aux;
}


double ELambda(double cutoff, double M1, double M2)
{
    double ELambdaAux = 0.5*( sqrt( pow(M1,2) + pow(cutoff,2) ) + sqrt( pow(M2,2) + pow(cutoff,2) ) );

    return ELambdaAux;
}


double pFunctionE(double M1, double M2, double E)
{   
    double sqrtArg = 16.0*pow(E,4) + pow(pow(M1,2)-pow(M2,2), 2) - 8.0*pow(E,2)*( pow(M1,2)+pow(M2,2) );

    double pFunctionEAux = sqrt( fabs(sqrtArg) )/( 4.0*E );

    if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
    { 
        cout << "Argument of the square root in pFunctionE is negative and larger then variable TFLI_ZERO!\n"; 
    }

    return pFunctionEAux;
}


double E1FunctionE(double M1, double M2, double E)
{   
    double sqrtArg = ( pow(-pow(M1,2) + pow(M2,2) + 4.0*pow(E,2), 2) )/pow(E,2);

    double E1FunctionEAux = 2.0*E - ( 1.0/4.0 )*sqrt( fabs(sqrtArg) );

    if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
    { 
        cout << "Argument of the square root in E1FunctionE is negative and larger then variable TFLI_ZERO!\n"; 
    }

    return E1FunctionEAux;
}


double E2FunctionE(double M1, double M2, double E)
{   
    double sqrtArg = ( pow(+pow(M1,2) - pow(M2,2) + 4.0*pow(E,2), 2) )/pow(E,2);

    double E2FunctionEAux = 2.0*E - ( 1.0/4.0 )*sqrt( fabs(sqrtArg) );

    if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
    { 
        cout << "Argument of the square root in E2FunctionE is negative and larger then variable TFLI_ZERO!\n"; 
    }

    return E2FunctionEAux;
}


double gPlusEtaZero3Momentum(double T, double effCP1, double effCP2, double M1, double M2, double eta, double E)
{
    double gPlusEtaAux = 1.0 - fermiDistribution(T, E1FunctionE(M1, M2, E) + eta*effCP1) - fermiDistribution(T, E2FunctionE(M1, M2, E) - eta*effCP2);

    return gPlusEtaAux;
}


double real16Pi2f1PairZero3Momentum3DCutoffNumerator(double E, void *parameters)
{   
    TwoFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double effCP1 = aux.getEffectiveChemicalPotential1();
    double effCP2 = aux.getEffectiveChemicalPotential2();
    double M1 = aux.getEffectiveMass1();
    double M2 = aux.getEffectiveMass2();
    double eta = aux.getEtaVariable();

    double real16Pi2f1PairNumerator = ( pFunctionE(M1, M2, E)/E )*gPlusEtaZero3Momentum(T, effCP1, effCP2, M1, M2, eta, E);

    return real16Pi2f1PairNumerator;
}


double real16Pi2f1PairZero3Momentum3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double integralPrecision)
{   
    int integrationWorkspace = 1000;

    //this is a temporary trick to perform some integrations, the case of equal masses and zero momentum has to be implemented
    if ( fabs(M1-M2)<TFLI_ZERO_MASS_DIFFERENCE ){ M2 = M2 + TFLI_ZERO; }

    double A = E0(M1, M2);
    double B = ELambda(cutoff, M1, M2);
    double sing = 0.5*( w + (-effCP2 + effCP1) );

    TwoFermionLine3DCutoffIntegrand aux(T, effCP1, effCP2, cutoff, M1, M2, w, 0.0);
    double real16Pi2f1Pairk0 = 0.0;

    ////////////////////////////////////////
    //eta = +1.0
    aux.setEtaVariable(+1.0);

    aux.setIntegralID("realf1PairPlusk0_A_B");
    Integration1DimGSLQAWCQAGS realf1PairPlusk0_A_B(A, B, -sing, &aux, real16Pi2f1PairZero3Momentum3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace);
    real16Pi2f1Pairk0 = real16Pi2f1Pairk0 + realf1PairPlusk0_A_B.evaluate();

    ////////////////////////////////////////
    //eta = -1.0
    aux.setEtaVariable(-1.0);

    aux.setIntegralID("realf1PairMinusk0_A_B");
    Integration1DimGSLQAWCQAGS realf1PairMinusk0_A_B(A, B, +sing, &aux, real16Pi2f1PairZero3Momentum3DCutoffNumerator, integralPrecision, integralPrecision, integrationWorkspace);
    real16Pi2f1Pairk0 = real16Pi2f1Pairk0 + realf1PairMinusk0_A_B.evaluate();

    return real16Pi2f1Pairk0;
}


double imag16Pi2f1PairZero3Momentum3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w)
{   
    double A = E0(M1, M2);
    double B = ELambda(cutoff, M1, M2);
    double sing = 0.5*( w + (-effCP2 + effCP1) );

    double imag16Pi2f1Pairk0 = 0.0;

    ////////////////////////////////////////
    //eta = +1.0
    double boxCarEtaPlus = heavisideTheta( -sing - A ) - heavisideTheta( -sing - B );
    if ( fabs(boxCarEtaPlus)>0 )
    {
        imag16Pi2f1Pairk0 = imag16Pi2f1Pairk0 + imag16Pi2f1PairSign(+1.0, symmetricImaginaryPart)*( M_PI*pFunctionE(M1, M2, -sing)/(-sing) )*gPlusEtaZero3Momentum(T, effCP1, effCP2, M1, M2, +1.0, -sing)*boxCarEtaPlus;
    }

    ////////////////////////////////////////
    //eta = -1.0
    double boxCarEtaMinus = heavisideTheta( sing - A ) - heavisideTheta( sing - B );
    if ( fabs(boxCarEtaMinus)>0 )
    {
        imag16Pi2f1Pairk0 = imag16Pi2f1Pairk0 + imag16Pi2f1PairSign(-1.0, symmetricImaginaryPart)*( M_PI*pFunctionE(M1, M2, +sing)/(+sing) )*gPlusEtaZero3Momentum(T, effCP1, effCP2, M1, M2, -1.0, +sing)*boxCarEtaMinus;
    }

    return imag16Pi2f1Pairk0;
}


////////////////////////////////////////////////////////////////////////////////////////
// calculate the special case in which the external three momentum is zero of the f1scat loop function
// lim k->0 f1scat: real16Pi2f1Scat3DCutoff and imag16Pi2f1Scat3DCutoff


double epsilon0(double M1, double M2)
{
    double epsilon0Aux = M2 - M1;

    return epsilon0Aux;
}


double epsilonLambda(double cutoff, double M1, double M2)
{
    double epsilonLambdaAux = sqrt( pow(M2,2) + pow(cutoff,2) ) - sqrt( pow(M1,2) + pow(cutoff,2) );

    return epsilonLambdaAux;
}


double pFunctionEpsilon(double M1, double M2, double epsilon)
{   
    double sqrtArg = (M1-M2-epsilon)*(M1+M2-epsilon)*(M1-M2+epsilon)*(M1+M2+epsilon)/( pow(2.0*epsilon,2) );

    double pFunctionEpsilonAux = sqrt( fabs(sqrtArg) );

    if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
    { 
        cout << "Argument of the square root in pFunctionEpsilon is negative and larger then variable TFLI_ZERO!\n"; 
    }

    return pFunctionEpsilonAux;
}


double E1FunctionEpsilon(double M1, double M2, double epsilon)
{   
    double sqrtArg = pow(-pow(M1,2) + pow(M2,2) + pow(epsilon,2),2)/pow(epsilon,2);

    double E1FunctionEpsilonAux = -epsilon + 0.5*sqrt( fabs(sqrtArg) );

    if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
    { 
        cout << "Argument of the square root in E1FunctionEpsilon is negative and larger then variable TFLI_ZERO!\n"; 
    }

    return E1FunctionEpsilonAux;
}


double E2FunctionEpsilon(double M1, double M2, double epsilon)
{   
    double sqrtArg = pow(+pow(M1,2) - pow(M2,2) + pow(epsilon,2),2)/pow(epsilon,2);

    double E2FunctionEpsilon = +epsilon + 0.5*sqrt( fabs(sqrtArg) );

    if ( sqrtArg<0 && fabs(sqrtArg)>TFLI_ZERO )
    { 
        cout << "Argument of the square root in E2FunctionEpsilon is negative and larger then variable TFLI_ZERO!\n"; 
    }

    return E2FunctionEpsilon;
}


double gMinusEtaZero3Momentum(double T, double effCP1, double effCP2, double M1, double M2, double eta, double epsilon)
{
    double gMinusEtaAux = + fermiDistribution(T, E1FunctionEpsilon(M1, M2, epsilon) - eta*effCP1) 
                          - fermiDistribution(T, E2FunctionEpsilon(M1, M2, epsilon) - eta*effCP2);

    return gMinusEtaAux;
}


double gMinusEtaZero3Momentum(double T, double effCP1, double effCP2, double eta, double E1, double E2)
{
    double gMinusEtaAux = fermiDistribution(T, E1 - eta*effCP1) - fermiDistribution(T, E2 - eta*effCP2);

    return gMinusEtaAux;
}


double real16Pi2f1ScatZero3Momentum3DCutoffDifferentMassesNumerator(double epsilon, void *parameters)
{   
    TwoFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double effCP1 = aux.getEffectiveChemicalPotential1();
    double effCP2 = aux.getEffectiveChemicalPotential2();
    double M1 = aux.getEffectiveMass1();
    double M2 = aux.getEffectiveMass2();
    double eta = aux.getEtaVariable();

    double real16Pi2f1ScatNumerator = ( 2.0*pFunctionEpsilon(M1, M2, epsilon)/epsilon )*gMinusEtaZero3Momentum(T, effCP1, effCP2, M1, M2, eta, epsilon);

    return real16Pi2f1ScatNumerator;
}


double real16Pi2f1ScatZero3Momentum3DCutoffDifferentMasses(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double integralPrecision)
{
    int integrationWorkspace = 1000;

    double a = epsilon0(M1, M2);
    double b = epsilonLambda(cutoff, M1, M2);
    double sing = ( w + (-effCP2 + effCP1) );

    TwoFermionLine3DCutoffIntegrand aux(T, effCP1, effCP2, cutoff, M1, M2, w, 0.0);
    double real16Pi2f1Scatk0 = 0.0;

    ////////////////////////////////////////
    //eta = +1.0
    aux.setEtaVariable(+1.0);

    aux.setIntegralID("realf1ScatPlusk0DiffM_a_b");
    Integration1DimGSLQAWCQAGS realf1ScatPlusk0DiffM_a_b(a, b, -sing, &aux, real16Pi2f1ScatZero3Momentum3DCutoffDifferentMassesNumerator, integralPrecision, integralPrecision, integrationWorkspace);
    real16Pi2f1Scatk0 = real16Pi2f1Scatk0 + realf1ScatPlusk0DiffM_a_b.evaluate();

    ////////////////////////////////////////
    //eta = -1.0
    aux.setEtaVariable(-1.0);

    aux.setIntegralID("realf1ScatMinusk0DiffM_a_b");
    Integration1DimGSLQAWCQAGS realf1ScatMinusk0DiffM_a_b(a, b, +sing, &aux, real16Pi2f1ScatZero3Momentum3DCutoffDifferentMassesNumerator, integralPrecision, integralPrecision, integrationWorkspace);
    real16Pi2f1Scatk0 = real16Pi2f1Scatk0 + realf1ScatMinusk0DiffM_a_b.evaluate();

    return real16Pi2f1Scatk0;
}


double real16Pi2f1ScatZero3Momentum3DCutoffEqualMassesIntegrand(double E, void *parameters)
{   
    TwoFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double effCP1 = aux.getEffectiveChemicalPotential1();
    double effCP2 = aux.getEffectiveChemicalPotential2();
    double M1 = aux.getEffectiveMass1();
    double M2 = aux.getEffectiveMass2();
    double w = aux.getZeroMomentum();

    if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
    {
        cout << "Masses should be equal in real16Pi2f1ScatZero3Momentum3DCutoffEqualMassesIntegrand! The difference is larger then TFLI_ZERO_MASS_DIFFERENCE!\n";
    }

    double sing = ( w + (-effCP2 + effCP1) );

    double real16Pi2f1ScatIntegrand = ( sqrt( pow(E,2) - pow(M1,2) )/E )*( + gMinusEtaZero3Momentum(T, effCP1, effCP2, +1.0, E, E)
                                                                           - gMinusEtaZero3Momentum(T, effCP1, effCP2, -1.0, E, E) )*( -2.0/sing );

    return real16Pi2f1ScatIntegrand;
}


double real16Pi2f1ScatZero3Momentum3DCutoffEqualMasses(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double integralPrecision)
{
    int integrationWorkspace = 1000;

    double a = M1;
    double b = sqrt( pow(cutoff,2) + pow(M1,2) );
    double sing = ( w + (-effCP2 + effCP1) );

    if ( fabs(sing)<=TFLI_ZERO ){ w = w + TFLI_ZERO; }

    TwoFermionLine3DCutoffIntegrand aux(T, effCP1, effCP2, cutoff, M1, M2, w, 0.0);
    double real16Pi2f1Scatk0 = 0.0;

    aux.setIntegralID("realf1ScatPlusk0EqualM_a_b");
    Integration1DimGSLQAGS realf1ScatPlusk0EqualM_a_b(a, b, &aux, real16Pi2f1ScatZero3Momentum3DCutoffEqualMassesIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
    real16Pi2f1Scatk0 = real16Pi2f1Scatk0 + realf1ScatPlusk0EqualM_a_b.evaluate();

    return real16Pi2f1Scatk0;
}


double imag16Pi2f1ScatZero3Momentum3DCutoffDifferentMasses(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w)
{   
    double a = epsilon0(M1, M2);
    double b = epsilonLambda(cutoff, M1, M2);
    double sing = ( w + (-effCP2 + effCP1) );

    double imag16Pi2f1Scatk0 = 0.0;

    ////////////////////////////////////////
    //eta = +1.0
    double boxCarEtaPlus = heavisideTheta( -sing - a ) - heavisideTheta( -sing - b );
    if ( fabs(boxCarEtaPlus)>0 )
    {
        imag16Pi2f1Scatk0 = imag16Pi2f1Scatk0 + imag16Pi2f1ScatSign(+1.0, -sing, symmetricImaginaryPart)*( 2.0*M_PI*pFunctionEpsilon(M1, M2, -sing)/(-sing) )*gMinusEtaZero3Momentum(T, effCP1, effCP2, M1, M2, +1.0, -sing)*boxCarEtaPlus;
    }

    ////////////////////////////////////////
    //eta = -1.0
    double boxCarEtaMinus = heavisideTheta( sing - a ) - heavisideTheta( sing - b );
    if ( fabs(boxCarEtaMinus)>0 )
    {
        imag16Pi2f1Scatk0 = imag16Pi2f1Scatk0 + imag16Pi2f1ScatSign(-1.0, +sing, symmetricImaginaryPart)*( 2.0*M_PI*pFunctionEpsilon(M1, M2, +sing)/(+sing) )*gMinusEtaZero3Momentum(T, effCP1, effCP2, M1, M2, -1.0, +sing)*boxCarEtaMinus;
    }

    return imag16Pi2f1Scatk0;
}


double real16Pi2f1Zero3Momentum3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double integralPrecision)
{   
    double real16Pi2f1 = real16Pi2f1PairZero3Momentum3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, integralPrecision);

    if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
    {
        real16Pi2f1 = real16Pi2f1 + real16Pi2f1ScatZero3Momentum3DCutoffDifferentMasses(T, effCP1, effCP2, cutoff, M1, M2, w, integralPrecision);
    }
    else
    {
        real16Pi2f1 = real16Pi2f1 + real16Pi2f1ScatZero3Momentum3DCutoffEqualMasses(T, effCP1, effCP2, cutoff, M1, M2, w, integralPrecision);
    }

    return real16Pi2f1;
}


double imag16Pi2f1Zero3Momentum3DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w)
{
    double imag16Pi2f1 = imag16Pi2f1PairZero3Momentum3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w);

    if ( fabs(M1-M2)>TFLI_ZERO_MASS_DIFFERENCE )
    {
        imag16Pi2f1 = imag16Pi2f1 + imag16Pi2f1ScatZero3Momentum3DCutoffDifferentMasses(T, effCP1, effCP2, cutoff, M1, M2, w);
    }
    else
    {
        imag16Pi2f1 = imag16Pi2f1 + 0.0;
    }

    return imag16Pi2f1;
}


double real16Pi2f13DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double k, double integralPrecision)
{   
    double real16Pi2f1 = 0.0;

    if ( k>=2*cutoff || k<0 )
    {   
        // zero by definition: k must be larger or equal to zero and if k>=2cutoff the spheres do not intersect
        real16Pi2f1 = 0.0;
    }
    else
    {
        if ( fabs(k)<1E-6 )
        {
            // k=0 case
            real16Pi2f1 = real16Pi2f1Zero3Momentum3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, integralPrecision);
        }
        else
        {
            // k finite case
            real16Pi2f1 = real16Pi2f1Finite3Momentum3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, k, integralPrecision);
        }
    }

    return real16Pi2f1;
}


double imag16Pi2f13DCutoff(double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double k)
{   
    double imag16Pi2f1 = 0.0;

    if ( k>=2*cutoff || k<0 )
    {   
        // zero by definition: k must be larger or equal to zero and if k>=2cutoff the spheres do not intersect
        imag16Pi2f1 = 0.0;
    }
    else
    {
        if ( fabs(k)<1E-6 )
        {
            // k=0 case
            imag16Pi2f1 = imag16Pi2f1Zero3Momentum3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w);
        }
        else
        {
            // k finite case
            imag16Pi2f1 = imag16Pi2f1Finite3Momentum3DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, k);
        }
    }

    return imag16Pi2f1;
}


gsl_complex klevanskyB0Integral3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double T, double effCP1, double effCP2, double cutoff, double M1, double M2, double w, double k, double integralPrecision)
{
    if ( reguScheme==cutoffOnDivergentIntegralsOnly )
    {
        cout << "The function klevanskyB0Integral3DCutoff is not defined for the NJL3DCutoffRegularizationScheme:cutoffOnDivergentIntegralsOnly! Aborting!\n";
        abort();
    }

    double ReB0Klev = real16Pi2f13DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, k, integralPrecision);
    double ImB0Klev = imag16Pi2f13DCutoff(T, effCP1, effCP2, cutoff, M1, M2, w, k);

    gsl_complex B0Klev = gsl_complex_rect(ReB0Klev, ImB0Klev);

    return B0Klev;
}
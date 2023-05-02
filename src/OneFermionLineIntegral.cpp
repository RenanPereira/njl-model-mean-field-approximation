#include <cmath>
#include <gsl/gsl_complex_math.h>
#include "OneFermionLineIntegral.h"
#include "TwoFermionLineIntegral.h"


/*
The relation between the f0 integral and the A integral in Klevansky notation is:
AKlevansky = - ( 16 Pi^2 )*f0
*/


//Primitive needed to calculate the f0 loop function
double f0Primitive(double p, double M)
{
    double f0Primitive = 0.5*p*sqrt( pow(p,2) + pow(M,2) ) - 0.5*pow(M,2)*log( p + sqrt( pow(p,2) + pow(M,2) ) );

    return f0Primitive;
}


double f0ConvergentIntegrand(double p, void *parameters)
{   
    OneFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double Cp = aux.getEffectiveChemicalPotential();
    double M = aux.getEffectiveMass();

    double E = sqrt( pow(p,2) + pow(M,2) );
    
    double f0Integrand = ( pow(p,2)/E )*( - fermiDistribution(T, E-Cp) - fermiDistribution(T, E+Cp) );

    return f0Integrand;
}


double f03DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double T, double Cp, double M, double integralPrecision)
{       
    int integrationWorkspace = 1000;
    OneFermionLine3DCutoffIntegrand params("f0ConvergentIntegral", T, Cp, M);

    double f0 = 0.0;
    if ( T>0.0 )
    {   
        //finite T
        f0 = ( 1.0/(4.0*M_PI*M_PI) )*( f0Primitive(cutoff, M) - f0Primitive(0, M) );

        if ( reguScheme==cutoffEverywhere || reguScheme==cutoffEverywhereWithCTmu )
        {   
            Integration1DimGSLQAGS f0Convergent(0.0, cutoff, &params, f0ConvergentIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
            f0 = f0 + ( 1.0/(4.0*M_PI*M_PI) )*f0Convergent.evaluate();
        }
        else if ( reguScheme==cutoffOnDivergentIntegralsOnly )
        {   
            Integration1DimGSLQAGIU f0Convergent(0.0, &params, f0ConvergentIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
            f0 = f0 + ( 1.0/(4.0*M_PI*M_PI) )*f0Convergent.evaluate();
        }
    }
    else
    {   
        //T=0 limit
        f0 = ( 1.0/(4.0*M_PI*M_PI) )*( f0Primitive(cutoff, M) - f0Primitive(fermiMomentum(Cp, M), M) );
    }

    return f0;
}


double realKlevanskyAIntegral3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double T, double Cp, double M, double integralPrecision)
{   
    double f0 = f03DCutoff(reguScheme, cutoff, T, Cp, M, integralPrecision);

    double AKlev = -(16.0*M_PI*M_PI)*f0;

    return AKlev;
}


gsl_complex klevanskyAIntegral3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double T, double Cp, double M, double integralPrecision)
{
    double ReAKlev = realKlevanskyAIntegral3DCutoff(reguScheme, cutoff, T, Cp, M, integralPrecision);
    double ImAKlev = 0.0;

    gsl_complex AKlev = gsl_complex_rect(ReAKlev, ImAKlev);

    return AKlev;
}


double sigmaNJL3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double Nc, double T, double Cp, double M, double integralPrecision)
{   
    /*
    //alternative way to calculate
    double f0 = f03DCutoff(reguScheme, cutoff, T, Cp, M, integralPrecision);
    double sigma = -4.0*Nc*M*f0;
    */
    
    double AKlev = realKlevanskyAIntegral3DCutoff(reguScheme, cutoff, T, Cp, M, integralPrecision);
    double sigma = ( Nc/(4.0*M_PI*M_PI) )*M*AKlev;

    return sigma;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//f0k loop function in the vacuum
double f0kVacuum(double cutoff, double M, double k)
{   
    double f0k = (1.0/2.0)*( + ( cutoff - k/2.0 )*sqrt( pow(cutoff - k/2.0, 2) + pow(M,2) )
                             + pow(M,2)*log( fabs(M) )
                             - pow(M,2)*log( cutoff - k/2.0 + sqrt( pow(cutoff - k/2.0, 2) + pow(M,2) ) ) 
                            );

    //k-dependent term, if k is to small, use series expansion around k~0 up to second order
    if (k>1E-10)
    {
        f0k = f0k + ( 2.0/(3.0*k) )*pow( pow(cutoff,2) - pow(k/2.0,2) + pow(M,2) , 3.0/2.0)
                  - ( 2.0/(3.0*k) )*( (cutoff-k/2.0)*(cutoff+k) + pow(M,2) )*sqrt( pow(cutoff - k/2.0, 2) + pow(M,2) );
    }
    else//series expansion
    {
        f0k = f0k + k*pow(cutoff,2)/( 4.0*sqrt( pow(cutoff,2) + pow(M,2) ) )
                  - k*k*( 3.0*pow(M,2)*cutoff + 2.0*pow(cutoff,3) )/( 12.0*pow( pow(cutoff,2) + pow(M,2) , 3.0/2.0) );
    }

    f0k = ( 1.0/(4.0*M_PI*M_PI) )*f0k;

    return f0k;
}


double f0kConvergentIntegrand1(double E, void *parameters)
{   
    OneFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double Cp = aux.getEffectiveChemicalPotential();
    double M = aux.getEffectiveMass();

    double f0kIntegrand = sqrt( pow(E,2) - pow(M,2) )*( - fermiDistribution(T, E-Cp) - fermiDistribution(T, E+Cp) );

    return f0kIntegrand;
}


double f0kConvergentIntegrand2(double E, void *parameters)
{   
    OneFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double Cp = aux.getEffectiveChemicalPotential();
    
    double f0kIntegrand = 1.0*( - fermiDistribution(T, E-Cp) - fermiDistribution(T, E+Cp) );

    return f0kIntegrand;
}


double f0kConvergentIntegrand3(double E, void *parameters)
{   
    OneFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double Cp = aux.getEffectiveChemicalPotential();
    
    double f0kIntegrand = ( -pow(E,2) )*( - fermiDistribution(T, E-Cp) - fermiDistribution(T, E+Cp) );

    return f0kIntegrand;
}


double f03DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double T, double Cp, double M, double k, double integralPrecision)
{   
    int integrationWorkspace = 1000;

    double f0k = 0.0;

    if ( k>=2.0*cutoff || k<0.0)
    {
        f0k = 0.0;
    }
    else
    {   
        //finite temperature
        //vacuum contribution (divergent)
        f0k = f0kVacuum(cutoff, M, k);

        //in-medium contribution (convergent)
        if ( reguScheme==cutoffEverywhere || reguScheme==cutoffEverywhereWithCTmu )
        {   
            //////////////////////////////////////////////////////////////////////
            //convergent term 1
            OneFermionLine3DCutoffIntegrand params1("f0kConvergent1", T, Cp, M);
            double min1 = M;
            double max1 = sqrt( pow(cutoff-k/2.0, 2) + pow(M,2) );
            
            Integration1DimGSLQAGS f0kConvergent1(min1, max1, &params1, f0kConvergentIntegrand1, integralPrecision, integralPrecision, integrationWorkspace);
            f0k = f0k + ( 1.0/(4.0*M_PI*M_PI) )*f0kConvergent1.evaluate();       

            if (k>1E-10)
            {   
                //////////////////////////////////////////////////////////////////////
                //convergent term 2
                OneFermionLine3DCutoffIntegrand params2("f0kConvergent2", T, Cp, M);
                double min2 = sqrt( pow(cutoff-k/2.0, 2) + pow(M,2) );
                double max2 = sqrt( pow(cutoff,2) - pow(k/2.0,2) + pow(M,2) );

                Integration1DimGSLQAGS f0kConvergent2(min2, max2, &params2, f0kConvergentIntegrand2, integralPrecision, integralPrecision, integrationWorkspace);
                f0k = f0k + ( 1.0/(4.0*M_PI*M_PI) )*( ( pow(cutoff,2) - pow(k/2.,2) + pow(M,2) )/k )*f0kConvergent2.evaluate();
                
                //////////////////////////////////////////////////////////////////////
                //convergent term 3
                OneFermionLine3DCutoffIntegrand params3("f0kConvergent3", T, Cp, M);
                double min3 = sqrt( pow(cutoff-k/2.0, 2) + pow(M,2) );
                double max3 = sqrt( pow(cutoff,2) - pow(k/2.0,2) + pow(M,2) );

                Integration1DimGSLQAGS f0kConvergent3(min3, max3, &params3, f0kConvergentIntegrand3, integralPrecision, integralPrecision, integrationWorkspace);
                f0k = f0k + ( 1.0/(4.0*M_PI*M_PI) )*( 1.0/k )*f0kConvergent3.evaluate();
            }     
        }
        else if ( reguScheme==cutoffOnDivergentIntegralsOnly )
        {   
            OneFermionLine3DCutoffIntegrand params1("f0kConvergent1", T, Cp, M);

            Integration1DimGSLQAGIU f0kConvergent1(M, &params1, f0kConvergentIntegrand1, integralPrecision, integralPrecision, integrationWorkspace);
            f0k = f0k + ( 1.0/(4.0*M_PI*M_PI) )*f0kConvergent1.evaluate();       
        }
    }

    return f0k;
}


double realKlevanskyAIntegral3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double T, double Cp, double M, double k, double integralPrecision)
{   
    double f0 = f03DCutoff(reguScheme, cutoff, T, Cp, M, k, integralPrecision);

    double AKlev = -(16.0*M_PI*M_PI)*f0;

    return AKlev;
}


gsl_complex klevanskyAIntegral3DCutoff(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double T, double Cp, double M, double k, double integralPrecision)
{
    double ReAKlev = realKlevanskyAIntegral3DCutoff(reguScheme, cutoff, T, Cp, M, k, integralPrecision);
    double ImAKlev = 0.0;

    gsl_complex AKlev = gsl_complex_rect(ReAKlev, ImAKlev);

    return AKlev;
}




////////////////////////////////////////////////////////////////////////////////////////


double gEEta(double T, double effCP1, double cutoff, double M1, double M2, double k, double eta, double E)
{   
    double C = dEdepsilon_epsilonMin(cutoff, M1, M2, k, E);
    double D = dEdepsilon_epsilonMax(cutoff, M1, M2, k, E);

/*
    double gPlusEtaAux;
    gEEtaAux =  + 0.5*( D - C )
                + 2.0*T*log( 1 + exp( (- E + C/2 - eta*effCP1 )/T ) )
                - 2.0*T*log( 1 + exp( (- E + D/2 - eta*effCP1 )/T ) );
*/

    //arguments inside the logarithms
    double argD1 = - E + 0.5*D - eta*effCP1;
    double argC1 = - E + 0.5*C - eta*effCP1;

    // + 0.5*( D - C )
    double gEEtaAux = 0.5*( D - C );
    
    // - 2.0*T*log( 1 + exp( (- E + D/2 - eta*effCP1 )/T ) )
    if ( isinf( exp(argD1/T) )==0 )
    { 
        gEEtaAux = gEEtaAux - 2.0*T*log( 1.0 + exp(argD1/T) ); 
    }
    else if( argD1>0 )
    { 
        gEEtaAux = gEEtaAux - 2.0*puiseuxExpansionTln1plusExpArgOverT(T, argD1); 
    }

    // + 2.0*T*log( 1 + exp( (- E + C/2 - eta*effCP1 )/T ) )
    if ( isinf( exp(argC1/T) )==0 )
    { 
        gEEtaAux = gEEtaAux + 2.0*T*log( 1.0 + exp(argC1/T) ); 
    }
    else if( argC1>0 )
    { 
        gEEtaAux = gEEtaAux + 2.0*puiseuxExpansionTln1plusExpArgOverT(T, argC1); 
    }

    return gEEtaAux;
}


double realKlevanskyAPair3DCutoffIntegrand(double E, void *parameters)
{   
    TwoFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double effCP1 = aux.getEffectiveChemicalPotential1();
    double cutoff = aux.getThreeMomentumCutoff();
    double M1 = aux.getEffectiveMass1();
    double M2 = aux.getEffectiveMass2();
    double k = aux.getThreeMomentum();
    double eta = aux.getEtaVariable();

    double realKlevanskyAPairIntegrand = E*gEEta(T, effCP1, cutoff, M1, M2, k, eta, E);

    return realKlevanskyAPairIntegrand;
}


double realKlevanskyAPair3DCutoff(double T, double effCP1, double cutoff, double M1, double M2, double k, double integralPrecision)
{   
    int integrationWorkspace = 1000;

    double A = dEdepsilon_EMin(cutoff, M1, M2, k);
    double B = dEdepsilon_EMax(cutoff, M1, M2, k);

    TwoFermionLine3DCutoffIntegrand aux(T, effCP1, 0.0, cutoff, M1, M2, 0.0, k);
    double realKlevanskyAPair = 0.0;


    ////////////////////////////////////////
    //eta = +1.0
    aux.setEtaVariable(+1.0);

    aux.setIntegralID("realKlevanskyAPairPlus");
    Integration1DimGSLQAGS realKlevanskyAPairPlus(A, B, &aux, realKlevanskyAPair3DCutoffIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
    realKlevanskyAPair = realKlevanskyAPair + realKlevanskyAPairPlus.evaluate();

    ////////////////////////////////////////
    //eta = -1.0
    aux.setEtaVariable(-1.0);

    aux.setIntegralID("realKlevanskyAPairMinus");
    Integration1DimGSLQAGS realKlevanskyAPairMinus(A, B, &aux, realKlevanskyAPair3DCutoffIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
    realKlevanskyAPair = realKlevanskyAPair + realKlevanskyAPairMinus.evaluate();

    realKlevanskyAPair = ( -2.0/k )*realKlevanskyAPair;

    return realKlevanskyAPair;
}


double gEpsilonEta(double T, double effCP1, double cutoff, double M1, double M2, double k, double eta, double epsilon)
{   
    double c = depsilondE_EMin(M1, M2, k, epsilon);
    double d = depsilondE_EMax(cutoff, M1, M2, k, epsilon);

/*
    double gPlusEtaAux;
    gEpsilonEtaAux =  + 0.5*( d - c )
                      - T*log( 1 + exp( ( epsilon/2 - c - eta*effCP1 )/T ) )
                      + T*log( 1 + exp( ( epsilon/2 - d - eta*effCP1 )/T ) );
*/

    //arguments inside the logarithms
    double argC1 = 0.5*epsilon - c - eta*effCP1;
    double argD1 = 0.5*epsilon - d - eta*effCP1;

    // + 0.5*( d - c )
    double gEpsilonEtaAux = 0.5*( d - c );
    
    // - T*log( 1 + exp( ( epsilon/2 - c - eta*effCP1 )/T ) )
    if ( isinf( exp(argC1/T) )==0 )
    { 
        gEpsilonEtaAux = gEpsilonEtaAux - T*log( 1.0 + exp(argC1/T) ); 
    }
    else if( argC1>0 )
    { 
        gEpsilonEtaAux = gEpsilonEtaAux - puiseuxExpansionTln1plusExpArgOverT(T, argC1); 
    }

    // + T*log( 1 + exp( ( epsilon/2 - d - eta*effCP1 )/T ) )
    if ( isinf( exp(argD1/T) )==0 )
    { 
        gEpsilonEtaAux = gEpsilonEtaAux + T*log( 1.0 + exp(argD1/T) ); 
    }
    else if( argD1>0 )
    { 
        gEpsilonEtaAux = gEpsilonEtaAux + puiseuxExpansionTln1plusExpArgOverT(T, argD1); 
    }

    return gEpsilonEtaAux;
}


double realKlevanskyAScat3DCutoffIntegrand(double epsilon, void *parameters)
{   
    TwoFermionLine3DCutoffIntegrand aux(parameters);
    double T = aux.getTemperature();
    double effCP1 = aux.getEffectiveChemicalPotential1();
    double cutoff = aux.getThreeMomentumCutoff();
    double M1 = aux.getEffectiveMass1();
    double M2 = aux.getEffectiveMass2();
    double k = aux.getThreeMomentum();
    double eta = aux.getEtaVariable();

    double realKlevanskyAScatIntegrand = epsilon*gEpsilonEta(T, effCP1, cutoff, M1, M2, k, eta, epsilon);

    return realKlevanskyAScatIntegrand;
}


double realKlevanskyAScat3DCutoff(double T, double effCP1, double cutoff, double M1, double M2, double k, double integralPrecision)
{   
    int integrationWorkspace = 1000;

    double a = depsilondE_epsilonMin(cutoff, M1, M2, k);
    double b = depsilondE_epsilonMax(cutoff, M1, M2, k);

    TwoFermionLine3DCutoffIntegrand aux(T, effCP1, 0.0, cutoff, M1, M2, 0.0, k);
    double realKlevanskyAScat = 0.0;


    ////////////////////////////////////////
    //eta = +1.0
    aux.setEtaVariable(+1.0);

    aux.setIntegralID("realKlevanskyAScatPlus");
    Integration1DimGSLQAGS realKlevanskyAScatPlus(a, b, &aux, realKlevanskyAScat3DCutoffIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
    realKlevanskyAScat = realKlevanskyAScat + realKlevanskyAScatPlus.evaluate();

    ////////////////////////////////////////
    //eta = -1.0
    aux.setEtaVariable(-1.0);

    aux.setIntegralID("realKlevanskyAScatMinus");
    Integration1DimGSLQAGS realKlevanskyAScatMinus(a, b, &aux, realKlevanskyAScat3DCutoffIntegrand, integralPrecision, integralPrecision, integrationWorkspace);
    realKlevanskyAScat = realKlevanskyAScat + realKlevanskyAScatMinus.evaluate();

    realKlevanskyAScat = ( -1.0/k )*realKlevanskyAScat;

    return realKlevanskyAScat;
}


double realKlevanskyA3DCutoffNEW(NJL3DCutoffRegularizationScheme reguScheme, double cutoff, double T, double Cp, double M, double k, double integralPrecision)
{   
    double realKlevanskyA = 0.0;

    if (k<0 || k>=2*cutoff)
    {   
        // zero by definition: k must be larger or equal to zero and if k>=2cutoff the spheres do not intersect
        realKlevanskyA = 0.0;
    }
    else
    {
        if ( fabs(k)<1E-6 )
        {
            // k=0 case
            realKlevanskyA = realKlevanskyAIntegral3DCutoff(reguScheme, cutoff, T, Cp, M, integralPrecision);
        }
        else
        {
            // k finite case
            
            //This quantity is independent of the value of M2 chosen!
            double M1 = M;
            double M2 = 1.1*M1;

            realKlevanskyA = realKlevanskyA + realKlevanskyAPair3DCutoff(T, Cp, cutoff, M1, M2, k, integralPrecision);
            realKlevanskyA = realKlevanskyA + realKlevanskyAScat3DCutoff(T, Cp, cutoff, M1, M2, k, integralPrecision);
        }
    }

    return realKlevanskyA;
}
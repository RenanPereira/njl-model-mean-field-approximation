#include <cmath>
#include <iostream>
#include "generalPhysicsAndMath.h"

using namespace std;


//Bose-Einstein distribution function
double boseDistribution(double temperature, double energy)
{
    double n = 1.0/( exp( energy/temperature ) - 1.0 );

    return n;
}


//Fermi-Dirac distribution function
double fermiDistribution(double temperature, double energy)
{
    double n = 1.0/( exp( energy/temperature ) + 1.0 );

    return n;
}


//Calculate Fermi's momentum given the mass and chemical potential
double fermiMomentum(double Cp, double M)
{   
    double lambdaF = 0.0;

    if ( pow(Cp,2) > pow(M,2) ){ lambdaF = sqrt( pow(Cp,2) - pow(M,2) ); } 
    else{ lambdaF = 0; }

    return lambdaF;
}


//Heaviside theta function, H(x). It is defined to be 0 if x<0 and 1 if x>=1
double heavisideTheta(double x)
{
    double heaviside = 0.0;
    if ( x<0 ){ heaviside = 0.0; }
    else{ heaviside = 1.0; }

    return heaviside;
}


//Puiseux series expansion of ln(1+x) at x->infinity: ln(1+x) ~ ln(x) + 1/x - 1/(2x^2) + 1/(3x^3) - 1/(4x^4) + O(1/x^5)
//In this case we implment: T ln( 1 + exp(A/T) ) ~ A + T*( exp(-A/T) - exp(-2*A/T)/2 + exp(-3*A/T)/3 - exp(-4*A/T)/4 )
double puiseuxExpansionTln1plusExpArgOverT(double T, double arg)
{
    double expansion = arg + T*( exp(-arg/T) - exp(-2.0*arg/T)/2.0 + exp(-3.0*arg/T)/3.0 - exp(-4.0*arg/T)/4.0 );

    return expansion;
}


//Function that returns the sign of x
double sign(double x)
{   
    if ( x>0 ){ return +1.0; }
    else
    {
        if ( x<0 ){ return -1.0; }
        else{ return 0.0; }
    }
}

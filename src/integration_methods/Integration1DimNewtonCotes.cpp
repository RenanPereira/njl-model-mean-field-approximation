#include <iostream>
#include <cmath>
#include "integration_methods/Integration1DimNewtonCotes.h"

using namespace std;


void Integration1DimNewtonCotes::setVariables(double lowerBoundAux, double upperBoundAux, int numberOfPartitionsAux, GeneralIntegrandParameters* integrandParametersAux, double integrandAux(double, void*), NewtonCotesRule ruleAux)
{
    lowerBound = lowerBoundAux;
    upperBound = upperBoundAux;
    numberOfPartitions = numberOfPartitionsAux;

    integrandParameters = integrandParametersAux;
    integrand = integrandAux;

    rule = ruleAux;

    //check if number of partitions are appropriate for the chosen rule
    if ( rule==trapezoidal )
    {
        if ( numberOfPartitions<2 )
        {
            cout << "Integration1DimNewtonCotes: to use the trapezoidal rule, at least 2 partitions are necessary!\n";
            abort();
        }
    }

    if ( rule==alternativeCompositeSimpson )
    {
        if ( numberOfPartitions<8 )
        {
            cout << "Integration1DimNewtonCotes: to use the alternative Composite Simpson rule, at least 8 partitions are necessary!\n";
            abort();
        }
    }
}


Integration1DimNewtonCotes::Integration1DimNewtonCotes(double lowerBoundAux, double upperBoundAux, int numberOfPartitionsAux, GeneralIntegrandParameters* integrandParametersAux, double integrandAux(double, void*))
{
    setVariables(lowerBoundAux, upperBoundAux, numberOfPartitionsAux, integrandParametersAux, integrandAux, trapezoidal);
}


Integration1DimNewtonCotes::Integration1DimNewtonCotes(double lowerBoundAux, double upperBoundAux, int numberOfPartitionsAux, GeneralIntegrandParameters* integrandParametersAux, double integrandAux(double, void*), NewtonCotesRule ruleAux)
{
    setVariables(lowerBoundAux, upperBoundAux, numberOfPartitionsAux, integrandParametersAux, integrandAux, ruleAux);
}


double Integration1DimNewtonCotes::evaluateTrapezoidal()
{   
    
    double dx = (upperBound-lowerBound)/(numberOfPartitions-1);
    double area = 0.0;
    for (int i = 1; i < numberOfPartitions-1; ++i)
    {   
        area = area + integrand(lowerBound + i*dx, integrandParameters);
    }

    area = dx*( area + 0.5*integrand(lowerBound, integrandParameters) + 0.5*integrand(upperBound, integrandParameters) );
    result = area;
    
    return area;
}


double Integration1DimNewtonCotes::evaluateAlternativeCompositeSimpson()
{   
    double dx = (upperBound-lowerBound)/(numberOfPartitions-1);
    double area = 0.0;
    for (int i = 4; i < numberOfPartitions-4; ++i)
    {   
        area = area + integrand(lowerBound + i*dx, integrandParameters);
    }

    area = ( dx/(48.0) )*( + 17.0*integrand(lowerBound + 0.0*dx, integrandParameters)
                           + 59.0*integrand(lowerBound + 1.0*dx, integrandParameters)
                           + 43.0*integrand(lowerBound + 2.0*dx, integrandParameters)
                           + 49.0*integrand(lowerBound + 3.0*dx, integrandParameters)
                           + 48.0*area
                           + 49.0*integrand(upperBound - 3.0*dx, integrandParameters)
                           + 43.0*integrand(upperBound - 2.0*dx, integrandParameters)
                           + 59.0*integrand(upperBound - 1.0*dx, integrandParameters)
                           + 17.0*integrand(upperBound - 0.0*dx, integrandParameters)
                         );
    result = area;

    return area;
}


double Integration1DimNewtonCotes::evaluate()
{   
    double area = 0.0;

    if ( rule==trapezoidal )
    {
        area = evaluateTrapezoidal();
    }
    else if ( rule==alternativeCompositeSimpson )
    {
        area = evaluateAlternativeCompositeSimpson();
    }

    return area;
}


double Integration1DimNewtonCotes::evaluateAvoidingSingularPoint(double singularity)
{ 

    double dx = (upperBound-lowerBound)/(numberOfPartitions-1);
    double delta = 2*dx;

    double area = 0.0;

    double a = 0;
    double b = 0;
    double orientation = +1;
    if ( lowerBound<upperBound )
    {
        a = lowerBound;
        b = upperBound;
        orientation = +1;
    }
    else
    {
        a = upperBound;
        b = lowerBound;
        orientation = -1;   
    }

    if ( fabs(b-singularity)<fabs(a-singularity) )
    {
        Integration1DimNewtonCotes LowSum(a, singularity-delta - fabs(b - (singularity+delta)), numberOfPartitions, integrandParameters, integrand, rule);
        area = area + LowSum.evaluate();
    
        Integration1DimNewtonCotes MiddleSum(singularity-delta - fabs(b - (singularity+delta)), singularity-delta, numberOfPartitions, integrandParameters, integrand, rule);
        area = area + MiddleSum.evaluate();

        Integration1DimNewtonCotes UpperSum(singularity+delta, b, numberOfPartitions, integrandParameters, integrand, rule);
        area = area + UpperSum.evaluate();

        area = orientation*area;
    }
    else
    {
        Integration1DimNewtonCotes LowSum(a, singularity-delta, numberOfPartitions, integrandParameters, integrand, rule);
        area = area + LowSum.evaluate();
        
        Integration1DimNewtonCotes MiddleSum(singularity+delta, (singularity+delta) + fabs(a - (singularity+delta)), numberOfPartitions, integrandParameters, integrand, rule);
        area = area + MiddleSum.evaluate();

        Integration1DimNewtonCotes UpperSum((singularity+delta) + fabs(a - (singularity+delta)), b, numberOfPartitions, integrandParameters, integrand, rule);
        area = area + UpperSum.evaluate();
    }
    
    return area;

}






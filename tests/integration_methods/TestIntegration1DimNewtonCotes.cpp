#include <iostream>
#include <cmath>

#include "TestIntegration1DimNewtonCotes.h"
#include "integration_methods/Integration1DimNewtonCotes.h"

double integrandTestNewtonCotes(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */
    double integrand = pow(x,2);

    return integrand;
}

double integrandTestNewtonCotesCPV(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */

    double integrand = ( 1.0 )/( x-1.0 );

    return integrand;
}

bool testIntegration1DimNewtonCotes(double relativeDifference)
{   
    cout << "Testing several Composite Trapezoidal Sum integration methods with different integrands.\n";
    cout << "All the integrals are normalized to 1.\n";

    double normalization = 0.0;

    TestIntegrandParameters aux1("integrandTestNewtonCotes");
    Integration1DimNewtonCotes newtonCotesSum(-1.0, +2.0, 200, &aux1, integrandTestNewtonCotes);
    normalization = (1.0/3.0);
    double resultNewtonCotesSum = normalization*newtonCotesSum.evaluate();
    cout << "resultNewtonCotesSum: " << resultNewtonCotesSum << "\n";

    TestIntegrandParameters aux2("integrandTestNewtonCotesCPV");
    Integration1DimNewtonCotes newtonCotesSumCPV(-1.0, 2.0, 100, &aux2, integrandTestNewtonCotesCPV, alternativeCompositeSimpson);
    normalization = (-1.0/log(2.0));
    double resultNewtonCotesSumCPV = normalization*newtonCotesSumCPV.evaluateAvoidingSingularPoint(1.0);
    cout << "resultNewtonCotesSumCPV: " << resultNewtonCotesSumCPV << "\n";

    bool testNewtonCotesSum = true;
    if ( fabs(resultNewtonCotesSum-1)>relativeDifference )
    {
        testNewtonCotesSum = false;
    }

    bool testNewtonCotesSumCPV = true;
    if ( fabs(resultNewtonCotesSumCPV-1)>relativeDifference )
    {
        testNewtonCotesSumCPV = false;
    }

    return testNewtonCotesSum && testNewtonCotesSumCPV;
}

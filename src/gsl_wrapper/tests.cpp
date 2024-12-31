#include <iostream>
#include <cmath>
#include "Integration1DimGSL.h"


double integrandTestGSL(double , void *);

double integrandTestGSLCauchy(double , void *);

double integrandTestGSLQAGP(double , void *);

double integrandTestGSLQAGI(double , void *);

double integrandTestGSLQAWS(double , void *);

bool hardcodedTestIntegration1DimGSL(double );


double integrandTestGSL(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */
    double integrand = pow(x,2);

    return integrand;
}


double integrandTestGSLCauchy(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */
    x = x; /* avoid unused parameter warning */

    double integrand = ( 1.0 );

    return integrand;
}


double integrandTestGSLQAGP(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */

    double integrand = (1.0/sqrt(x))*log(x);

    return integrand;
}


double integrandTestGSLQAGI(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */

    double integrand = exp( -pow(x,2) );

    return integrand;
}


double integrandTestGSLQAWS(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */

    double integrand = pow(x,2);

    return integrand;
}


bool hardcodedTestIntegration1DimGSL(double relativeDifference)
{   
    cout << "Testing several GSL integration methods with different integrands.\n";
    cout << "All the integrals are normalized to 1.\n";

    double normalization = 0.0;

    TestIntegrandParameters aux1("integrandTestGSL");
    Integration1DimGSLQNG integralQNG(-1.0, +2.0, &aux1, integrandTestGSL, 1E-8, 1E-8);
    normalization = (1.0/3.0);
    double resultQNG = normalization*integralQNG.evaluate();
    cout << "resultQNG: " << resultQNG << "\n";


    Integration1DimGSLQAG integralQAG(-1.0, +2.0, &aux1, integrandTestGSL, 1E-8, 1E-8, 1000, 1);
    normalization = (1.0/3.0);
    double resultQAG = normalization*integralQAG.evaluate();
    cout << "resultQAG: " << resultQAG << "\n";


    Integration1DimGSLQAGS integralQAGS(-1.0, +2.0, &aux1, integrandTestGSL, 1E-8, 1E-8, 1000);
    normalization = (1.0/3.0);
    double resultQAGS = normalization*integralQAGS.evaluate();
    cout << "resultQAGS: " << resultQAGS << "\n";


    TestIntegrandParameters aux2("integrandTestGSLCauchy");
    Integration1DimGSLQAWC integralQAWC(-1.0, +2.0, 1.0, &aux2, integrandTestGSLCauchy, 1E-8, 1E-8, 1000);
    normalization = (-1.0/log(2.0));
    double resultQAWC = normalization*integralQAWC.evaluate();
    cout << "resultQAWC: " << resultQAWC << "\n";


    TestIntegrandParameters aux3("integrandTestGSLQAGP");
    Integration1DimGSLQAGP integralQAGP(0, +1.0, {0.0}, &aux3, integrandTestGSLQAGP, 1E-8, 1E-8, 1000);
    normalization = -0.25;
    double resultQAGP = normalization*integralQAGP.evaluate();
    cout << "resultQAGP: " << resultQAGP << "\n";

    Integration1DimGSLCQUAD integralCQUAD(0, +1.0, &aux3, integrandTestGSLQAGP, 1E-8, 1E-8, 1000);
    normalization = -0.25;
    double resultCQUAD = normalization*integralCQUAD.evaluate();
    cout << "resultCQUAD: " << resultCQUAD << "\n";


    TestIntegrandParameters aux4("integrandTestGSLQAGI");
    Integration1DimGSLQAGI integralQAGI(&aux4, integrandTestGSLQAGI, 1E-8, 1E-8, 1000);
    normalization = ( sqrt(1.0/M_PI) );
    double resultQAGI = normalization*integralQAGI.evaluate();
    cout << "resultQAGI: " << resultQAGI << "\n";

    Integration1DimGSLQAGIU integralQAGIU(0, &aux4, integrandTestGSLQAGI, 1E-8, 1E-8, 1000);
    normalization = ( sqrt(1.0/M_PI) );
    double resultQAGIU = normalization*2*integralQAGIU.evaluate();
    cout << "resultQAGIU: " << resultQAGIU << "\n";

    Integration1DimGSLQAGIL integralQAGIL(0, &aux4, integrandTestGSLQAGI, 1E-8, 1E-8, 1000);
    normalization = ( sqrt(1.0/M_PI) );
    double resultQAGIL = normalization*2*integralQAGIL.evaluate();
    cout << "resultQAGIL: " << resultQAGIL << "\n";

    TestIntegrandParameters aux5("integrandTestGSLQAWS");
    Integration1DimGSLQAWS integralQAWS(-1.0, +2.0, &aux1, integrandTestGSLQAWS, 1E-8, 1E-8, 1000, -0.5, 0.0, 0, 0);
    normalization = ( 5.0/( 8.0*sqrt(3) ) );
    double resultQAWS = normalization*integralQAWS.evaluate();
    cout << "resultQAWS: " << resultQAWS << "\n";


    TestIntegrandParameters aux6("integrandTestGSLQAWCQAGS");
    Integration1DimGSLQAWCQAGS integralQAWCQAGSIn(-5.0, 5.0, 4.0, &aux6, integrandTestGSLCauchy, 1E-8, 1E-8, 1000);
    normalization = (-1.0/log(9.0));
    
    double resultQAWCQAGSIn = normalization*integralQAWCQAGSIn.evaluate();
    cout << "resultQAWCQAGSIn: " << resultQAWCQAGSIn << "\n";

    double compositeSumQAWCQAGSIn = normalization*integralQAWCQAGSIn.evaluateIntegration1DimNewtonCotes(10, alternativeCompositeSimpson);
    cout << "compositeSumQAWCQAGSIn: " << compositeSumQAWCQAGSIn << "\n";

    Integration1DimGSLQAWCQAGS integralQAWCQAGSOut(-5.0, 5.0, 8.0, &aux6, integrandTestGSLCauchy, 1E-8, 1E-8, 1000);
    normalization = (-1.0/log(13.0/3.0));
    
    double resultQAWCQAGSOut = normalization*integralQAWCQAGSOut.evaluate();
    cout << "resultQAWCQAGSOut: " << resultQAWCQAGSOut << "\n";

    double newtonCotesSumQAWCQAGSOut = normalization*integralQAWCQAGSOut.evaluateIntegration1DimNewtonCotes(10, alternativeCompositeSimpson);
    cout << "newtonCotesSumQAWCQAGSOut: " << newtonCotesSumQAWCQAGSOut << "\n";


    bool testIntegralQNG = true;
    if ( fabs(resultQNG-1)>relativeDifference )
    { 
        testIntegralQNG = false; 
    }

    bool testIntegralQAG = true;
    if ( fabs(resultQAG-1)>relativeDifference )
    { 
        testIntegralQAG = false; 
    }

    bool testIntegralQAGS = true;
    if ( fabs(resultQAGS-1)>relativeDifference )
    { 
        testIntegralQAGS = false; 
    }

    bool testIntegralQAWC = true;
    if ( fabs(resultQAWC-1)>relativeDifference )
    { 
        testIntegralQAWC = false; 
    }

    bool testIntegralQAGP = true;
    if ( fabs(resultQAGP-1)>relativeDifference )
    { 
        testIntegralQAGP = false; 
    }

    bool testIntegralCQUAD = true;
    if ( fabs(resultCQUAD-1)>relativeDifference )
    { 
        testIntegralCQUAD = false; 
    }

    bool testIntegralQAGI = true;
    if ( fabs(resultQAGI-1)>relativeDifference )
    { 
        testIntegralQAGI = false; 
    }

    bool testIntegralQAGIU = true;
    if ( fabs(resultQAGIU-1)>relativeDifference )
    { 
        testIntegralQAGIU = false; 
    }

    bool testIntegralQAGIL = true;
    if ( fabs(resultQAGIL-1)>relativeDifference )
    { 
        testIntegralQAGIL = false; 
    }

    bool testIntegralQAWS = true;
    if ( fabs(resultQAWS-1)>relativeDifference )
    { 
        testIntegralQAWS = false; 
    }

    bool testIntegralQAWCQAGSIn = true;
    if ( fabs(resultQAWCQAGSIn-1)>relativeDifference )
    { 
        testIntegralQAWCQAGSIn = false; 
    }

    bool testIntegralQAWCQAGSOut = true;
    if ( fabs(resultQAWCQAGSOut-1)>relativeDifference )
    { 
        testIntegralQAWCQAGSOut = false; 
    }

    return testIntegralQNG &&
           testIntegralQAG &&
           testIntegralQAGS &&
           testIntegralQAWC &&
           testIntegralQAGP &&
           testIntegralCQUAD &&
           testIntegralQAGI &&
           testIntegralQAGIU &&
           testIntegralQAGIL &&
           testIntegralQAWS &&
           testIntegralQAWCQAGSIn &&
           testIntegralQAWCQAGSOut;
}


int main() 
{
    bool allTestsPassed = true;

    // Run each test and collect the result
    allTestsPassed &= hardcodedTestIntegration1DimGSL(1E-8);

    if (allTestsPassed) 
    {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } 
    else 
    {
        std::cout << "Some tests failed!" << std::endl;
        return 1;
    }
}

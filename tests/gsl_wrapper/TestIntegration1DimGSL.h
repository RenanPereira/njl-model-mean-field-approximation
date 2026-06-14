#ifndef TESTINTEGRATION1DIMGSL_H
#define TESTINTEGRATION1DIMGSL_H

class TestIntegration1DimGSL
{
private:
    static double integrandTestGSL(double , void *);
    static double integrandTestGSLCauchy(double , void *);
    static double integrandTestGSLQAGP(double , void *);
    static double integrandTestGSLQAGI(double , void *);
    static double integrandTestGSLQAWS(double , void *);

public:
    static bool hardcodedTestIntegration1DimGSL(double );
};

#endif

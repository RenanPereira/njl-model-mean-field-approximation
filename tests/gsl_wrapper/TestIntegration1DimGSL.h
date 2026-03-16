#ifndef TESTINTEGRATION1DIMGSL_H
#define TESTINTEGRATION1DIMGSL_H

double integrandTestGSL(double , void *);

double integrandTestGSLCauchy(double , void *);

double integrandTestGSLQAGP(double , void *);

double integrandTestGSLQAGI(double , void *);

double integrandTestGSLQAWS(double , void *);

bool hardcodedTestIntegration1DimGSL(double );

#endif

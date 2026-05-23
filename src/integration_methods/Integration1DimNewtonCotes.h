#ifndef INTEGRATION1DIMNEWTONCOTES_H
#define INTEGRATION1DIMNEWTONCOTES_H

#include <iostream>
#include "integration_methods/GeneralIntegrandParameters.h"


enum class NewtonCotesRule
{ 
	TRAPEZOIDAL,
	ALTERNATIVE_COMPOSITE_SIMPSON
};

class Integration1DimNewtonCotes
{	
public:
	double lowerBound;
	double upperBound;
	int numberOfPartitions;

	GeneralIntegrandParameters* integrandParameters;
	double (* integrand)(double, void * parameters);

	NewtonCotesRule rule;
	
	double result;

public:
	Integration1DimNewtonCotes(){};
	Integration1DimNewtonCotes(double , double , int , GeneralIntegrandParameters* , double (double, void*) );
	Integration1DimNewtonCotes(double , double , int , GeneralIntegrandParameters* , double (double, void*) , NewtonCotesRule );

	void setVariables(double , double , int , GeneralIntegrandParameters* , double (double, void*), NewtonCotesRule );

	double evaluateTrapezoidal();
	double evaluateAlternativeCompositeSimpson();
	double evaluate();
	double evaluateAvoidingSingularPoint(double singularity);
};


class TestIntegrandParameters : public GeneralIntegrandParameters
{
public:
    std::string integralID;

public:
    TestIntegrandParameters(){};
    TestIntegrandParameters(std::string integralIDAux){ integralID = integralIDAux; };

    void printIntegrandVariables() override
    {   
        std::cout << integralID << "\n";
    }
    
};


#endif

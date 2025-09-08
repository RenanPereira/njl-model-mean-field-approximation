#ifndef INTEGRATION1DIMNEWTONCOTES_H
#define INTEGRATION1DIMNEWTONCOTES_H

#include <iostream>

using namespace std;


enum NewtonCotesRule { trapezoidal, alternativeCompositeSimpson };


class GeneralIntegrandParameters
{
public:
    virtual void printIntegrandVariables()
    {
        cout << "The method 'printIntegrandVariables' is using the default method from the class GeneralIntegrandParameters!" << "\n";
    }

    virtual void behaviourAfterFailedIntegration()
    { 
    	cout << "Using the default behaviour after a failed integration: abort!" << "\n";
    	abort(); 
    };
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
    string integralID;

public:
    TestIntegrandParameters(){};
    TestIntegrandParameters(string integralIDAux){ integralID = integralIDAux; };

    void printIntegrandVariables() override
    {   
        cout << integralID << "\n";
    }
    
};


#endif

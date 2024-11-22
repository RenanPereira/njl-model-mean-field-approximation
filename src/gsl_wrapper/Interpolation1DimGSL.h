#ifndef INTERPOLATION1DIMGSL_H
#define INTERPOLATION1DIMGSL_H

#include <vector> 
#include <iostream>
#include <gsl/gsl_spline.h>
#include "OneVariableFunction.h"
#include "root_solver_gsl.h"

using namespace std;


enum InterpolationMethod { linear, steffen, cubic, akima, polynomial };

class Interpolation1DimGSL
{
private:
	InterpolationMethod method;
    vector<double> discretizedVariable;
    vector<double> discretizedFunction;
    int interpolationSize;
    double interpolationLowerBound;
    double interpolationUpperBound;
    gsl_spline *spline;
    gsl_interp_accel *accelerator;

private:
    void setSpline();
    void unsetSpline();
    void setAccelerator();
    void unsetAccelerator();
    
public:
    Interpolation1DimGSL(InterpolationMethod , vector<double> , vector<double> );//constructor
    Interpolation1DimGSL(InterpolationMethod, double , double , int , OneVariableFunction );
    
    ~Interpolation1DimGSL()//destructor
    {   
        unsetSpline();
        unsetAccelerator();
    }

    Interpolation1DimGSL(const Interpolation1DimGSL &interpolation)//copy constructor
    {   
        method = interpolation.method;
        discretizedVariable = interpolation.discretizedVariable;
        discretizedFunction = interpolation.discretizedFunction;
        interpolationSize = interpolation.interpolationSize;
        interpolationLowerBound = interpolation.interpolationLowerBound;
        interpolationUpperBound = interpolation.interpolationUpperBound;
        setSpline();
        setAccelerator();
    }    

    Interpolation1DimGSL& operator=(const Interpolation1DimGSL& interpolation)//copy assignment
    {   
        if (this == &interpolation){ return *this; }
 
        method = interpolation.method;
        discretizedVariable = interpolation.discretizedVariable;
        discretizedFunction = interpolation.discretizedFunction;
        interpolationSize = interpolation.interpolationSize;
        interpolationLowerBound = interpolation.interpolationLowerBound;
        interpolationUpperBound = interpolation.interpolationUpperBound;
        unsetSpline();
        unsetAccelerator();
        setSpline();
        setAccelerator();

        return *this;
    }

    double getInterpolationLowerBound(){ return interpolationLowerBound; }
    double getInterpolationUpperBound(){ return interpolationUpperBound; }
    bool valueInInterpolationInterval(double );

    double evaluate(double );
    double evaluate1stDerivative(double );
    double evaluate2ndDerivative(double );
    double evaluateIntegral(double , double );
    
    vector<double> findRoots(RootFindingMethod, double );
    vector<double> findRoots1stDerivative(RootFindingMethod , double );
    vector<double> findRoots2ndDerivative(RootFindingMethod , double );
    
    void tests(OneVariableFunction );
};

double polynomialTest(double , void* );

double equationToFindRootsOfInterpolation(double , void *);

double equationToFindRootsOfInterpolation1stDerivative(double , void *);

double equationToFindRootsOfInterpolation2ndDerivative(double , void *);


#endif
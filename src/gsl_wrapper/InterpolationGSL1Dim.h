#ifndef INTERPOLATIONGSL1DIM_H
#define INTERPOLATIONGSL1DIM_H

#include <vector> 
#include <gsl/gsl_spline.h>
#include "math_utils/OneVariableFunction.h"
#include "gsl_wrapper/root_solver_gsl.h"

enum InterpolationGSL1DimMethod { linear, steffen, cubic, akima, polynomial };

class InterpolationGSL1Dim
{
private:
	InterpolationGSL1DimMethod method;
    std::vector<double> discretizedVariable;
    std::vector<double> discretizedFunction;
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
    InterpolationGSL1Dim(InterpolationGSL1DimMethod , std::vector<double> , std::vector<double> );//constructor
    InterpolationGSL1Dim(InterpolationGSL1DimMethod, double , double , int , OneVariableFunction );
    
    ~InterpolationGSL1Dim()//destructor
    {   
        unsetSpline();
        unsetAccelerator();
    }

    InterpolationGSL1Dim(const InterpolationGSL1Dim &interpolation)//copy constructor
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

    InterpolationGSL1Dim& operator=(const InterpolationGSL1Dim& interpolation)//copy assignment
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
    
    std::vector<double> findRoots(RootFindingMethod, double );
    std::vector<double> findRoots1stDerivative(RootFindingMethod , double );
    std::vector<double> findRoots2ndDerivative(RootFindingMethod , double );
    
    void tests(OneVariableFunction );
};

double polynomialTest(double , void* );

double equationToFindRootsOfInterpolation(double , void *);

double equationToFindRootsOfInterpolation1stDerivative(double , void *);

double equationToFindRootsOfInterpolation2ndDerivative(double , void *);


#endif
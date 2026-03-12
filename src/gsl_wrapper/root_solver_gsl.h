#ifndef ROOT_SOLVER_GSL_H
#define ROOT_SOLVER_GSL_H

#include <vector>
#include <map>
#include <gsl/gsl_vector.h>
#include "gsl_wrapper/ComplexSquareMatrixGSL.h"


// Maximum number of iterations in the root-finding loops
inline constexpr int MAX_ITERATIONS = 1000;

enum MultiRootFindingMethod 
{ 
    HYBRIDS, 
    HYBRID, 
    DNEWTON, 
    BROYDEN
};

inline const std::map<MultiRootFindingMethod, std::string> MultiRootFindingMethodMap = 
{
    {MultiRootFindingMethod::HYBRIDS, "HYBRIDS"},
    {MultiRootFindingMethod::HYBRID, "HYBRID"},
    {MultiRootFindingMethod::DNEWTON, "DNEWTON"},
    {MultiRootFindingMethod::BROYDEN, "BROYDEN"}
};

enum RootFindingMethod 
{ 
    brent, 
    bisection, 
    falsepos 
};

std::string toString(MultiRootFindingMethod );

MultiRootFindingMethod stringToMultiRootFindingMethod(const std::string& );

bool isValidMultiRootFindingMethod(const std::string& );

void multiDimensionalRootFind(int , double , double* , void* , int placeholder_f(const gsl_vector*, void*, gsl_vector*), MultiRootFindingMethod );

double OneDimensionalRootFind(double , double , double , void* , double placeholder_f(double, void*), RootFindingMethod );

std::vector<double> multiDimensionalRootFindRelativeErrors(int , double* , void* , int placeholder_f(const gsl_vector*, void*, gsl_vector*));

int multiDimensionalRootFindTestResidual(int , double , double* , void* , int placeholder_f(const gsl_vector*, void*, gsl_vector*));

std::vector<gsl_complex> sortGSLComplexNumbersByAbsoluteSize(std::vector<gsl_complex> );

gsl_complex cardanoA(gsl_complex , gsl_complex , gsl_complex , gsl_complex );

gsl_complex cardanoB(gsl_complex , gsl_complex , gsl_complex );

std::vector<gsl_complex> solveCubicEquationCardano(gsl_complex , gsl_complex , gsl_complex , gsl_complex );

std::vector<gsl_complex> calculateEigenvalues3By3ComplexMatrix(ComplexSquareMatrixGSL );

#endif
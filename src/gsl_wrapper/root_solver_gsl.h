#ifndef ROOT_SOLVER_GSL_H
#define ROOT_SOLVER_GSL_H

#include <vector>
#include <map>
#include <gsl/gsl_vector.h>
#include "gsl_wrapper/ComplexSquareMatrixGSL.h"


// Maximum number of iterations in the root-finding loops
const int MAX_ITERATIONS = 1000;

enum MultiRootFindingMethod 
{ 
    hybrids, 
    hybrid, 
    dnewton, 
    broyden
};

static const std::map<MultiRootFindingMethod, std::string> MultiRootFindingMethodMap = 
{
    {MultiRootFindingMethod::hybrids, "hybrids"},
    {MultiRootFindingMethod::hybrid, "hybrid"},
    {MultiRootFindingMethod::dnewton, "dnewton"},
    {MultiRootFindingMethod::broyden, "broyden"}
};

string toStringMultiRootFindingMethod(MultiRootFindingMethod );

MultiRootFindingMethod stringToMultiRootFindingMethod(const string& );

bool isValidMultiRootFindingMethod(const string& );

void multiDimensionalRootFind(int , double , double* , void* , int (const gsl_vector*, void*, gsl_vector*), MultiRootFindingMethod );

enum RootFindingMethod { brent, bisection, falsepos };

double OneDimensionalRootFind(double , double , double , void* , double placeholder_f (double, void*), RootFindingMethod );

vector<double> multiDimensionalRootFindRelativeErrors(int , double* , void* , int (const gsl_vector*, void*, gsl_vector*));

int multiDimensionalRootFindTestResidual(int , double , double* , void* , int (const gsl_vector*, void*, gsl_vector*));

vector<gsl_complex> sortGSLComplexNumbersByAbsoluteSize(vector<gsl_complex> );

gsl_complex cardanoA(gsl_complex , gsl_complex , gsl_complex , gsl_complex );

gsl_complex cardanoB(gsl_complex , gsl_complex , gsl_complex );

vector<gsl_complex> solveCubicEquationCardano(gsl_complex , gsl_complex , gsl_complex , gsl_complex );

vector<gsl_complex> calculateEigenvalues3By3ComplexMatrix(ComplexSquareMatrixGSL );

#endif
#ifndef ROOT_SOLVER_GSL_H
#define ROOT_SOLVER_GSL_H

#include <vector>
#include <gsl/gsl_vector.h>
#include "ComplexSquareMatrixGSL.h"

using namespace std;

enum MultiRootFindingMethod { hybrids, hybrid, dnewton, broyden };

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

double linearFit(double , double , double , double , double );

#endif
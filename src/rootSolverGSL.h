#ifndef ROOTSOLVERGSL_H
#define ROOTSOLVERGSL_H

#include <vector>
#include <gsl/gsl_vector.h>

using namespace std;

enum MultiRootFindingMethod { hybrids, hybrid, dnewton, broyden };

void multiDimensionalRootFind(int , double , double* , void* , int (const gsl_vector*, void*, gsl_vector*), MultiRootFindingMethod );

enum RootFindingMethod { brent, bisection, falsepos };

double OneDimensionalRootFind(double , double , double , void* , double placeholder_f (double, void*), RootFindingMethod );

vector<double> multiDimensionalRootFindRelativeErrors(int , double* , void* , int (const gsl_vector*, void*, gsl_vector*));

int multiDimensionalRootFindTestResidual(int , double , double* , void* , int (const gsl_vector*, void*, gsl_vector*));

#endif
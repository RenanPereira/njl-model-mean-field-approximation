#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <iostream>
#include "rootSolverGSL.h"

using namespace std;


//Multi-dimensional root-finding
void multiDimensionalRootFind(int n_eqs, double precision, double* x_init, void* params, int placeholder_f(const gsl_vector*, void*, gsl_vector*), MultiRootFindingMethod method)
{
	//set number of equations
	const size_t dim = n_eqs;

	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;
	
	//choose root fiding method
	if ( method==hybrids )
	{
		T = gsl_multiroot_fsolver_hybrids;
	}
	else if ( method==hybrid )
	{
		T = gsl_multiroot_fsolver_hybrid;
	}
	else if ( method==dnewton )
	{
		T = gsl_multiroot_fsolver_dnewton;
	}
	else if ( method==broyden )
	{
		T = gsl_multiroot_fsolver_broyden;
	}

	//allocate memory
	s = gsl_multiroot_fsolver_alloc(T,dim);		

	gsl_vector *x = gsl_vector_alloc(dim);

	//set initial guess
	for (int i = 0; i < n_eqs; ++i) gsl_vector_set(x,i,x_init[i]);	
	
	//set system of equations that will enter the method
	gsl_multiroot_function f = {placeholder_f, dim, params};

	//set up method
	gsl_multiroot_fsolver_set (s, &f, x);

	//start iteration
	int counter = 0;
	int status = GSL_CONTINUE;

	while (status == GSL_CONTINUE && counter < 1000)
	{
		++counter;
		status = gsl_multiroot_fsolver_iterate(s);

		if (status) break;

		status = gsl_multiroot_test_residual (s->f, precision);
	}

	//set solution to initial vector
	for (int i = 0; i < n_eqs; ++i) x_init[i] = gsl_vector_get(s->x, i);

	//free memory
	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
}


//One-dimensional root-finding
double OneDimensionalRootFind(double precision, double x_low, double x_high, void* params, double placeholder_f (double, void*), RootFindingMethod method)
{	
	gsl_function F;
    F.function = placeholder_f;
    F.params = params;

	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	
	//choose root fiding method
	if ( method==brent )
	{
		T = gsl_root_fsolver_brent;
	}
	else if ( method==bisection )
	{
		T = gsl_root_fsolver_bisection;
	}
	else if ( method==falsepos )
	{
		T = gsl_root_fsolver_falsepos;
	}

	//allocate memory
	s =  gsl_root_fsolver_alloc (T);		

	gsl_root_fsolver_set (s, &F, x_low, x_high);

	//start iteration
	int counter = 0;
	int status = GSL_CONTINUE;

	double root;
	while (status == GSL_CONTINUE && counter < 1000){

		++counter;
		status = gsl_root_fsolver_iterate(s);

		root = gsl_root_fsolver_root (s);
		x_low = gsl_root_fsolver_x_lower (s);
		x_high = gsl_root_fsolver_x_upper (s);

		if (status) break;

		status = gsl_root_test_interval (x_low, x_high, 0, precision);
	}

	//free memory
	gsl_root_fsolver_free (s);

	return root;
}


//Return the relative errors given a set of "roots" (x) and the multi root system of equations placeholder_f
vector<double> multiDimensionalRootFindRelativeErrors(int n, double* x, void* params, int placeholder_f(const gsl_vector*, void*, gsl_vector*))
{
    //create gsl vectors
    gsl_vector *xAux = gsl_vector_alloc(n);
    gsl_vector *fAux = gsl_vector_alloc(n);
    for (int i = 0; i < n; ++i)
    {
        gsl_vector_set(xAux, i, x[i]);
        gsl_vector_set(fAux, i, 0.00);//fAux is set to zero by default
    }

    //use the values of the x array in the system of equations placeholder_f, the result is stored in fAux
    placeholder_f(xAux, params, fAux);

    //create a vector to return the results
    vector<double> relativeErrors(n);
    for (int i = 0; i < n; ++i){ relativeErrors[i] = gsl_vector_get(fAux,i); }

    return relativeErrors;
}


int multiDimensionalRootFindTestResidual(int n, double precision, double* x, void* params, int placeholder_f(const gsl_vector*, void*, gsl_vector*))
{
    //create gsl vectors
    gsl_vector *xAux = gsl_vector_alloc(n);
    gsl_vector *fAux = gsl_vector_alloc(n);
    for (int i = 0; i < n; ++i)
    {
        gsl_vector_set(xAux, i, x[i]);
        gsl_vector_set(fAux, i, 0.00);//fAux is set to zero by default
    }

    //use the values of the x array in the system of equations placeholder_f, the result is stored in fAux
    placeholder_f(xAux, params, fAux);

    return gsl_multiroot_test_residual(fAux, precision);
}

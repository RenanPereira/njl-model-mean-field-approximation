#include <iostream>
#include <tuple>
#include <algorithm>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "gsl_wrapper/root_solver_gsl.h"

using namespace std;


string toString(MultiRootFindingMethod method) 
{
    // Check if the method exists in the map using count
    if (MultiRootFindingMethodMap.count(method))
    {
        return MultiRootFindingMethodMap.at(method);
    } 
    else 
    {
        cout << "Error: MultiRootFindingMethod not found in map! Returning UNKNOWN." << endl;
        return "UNKNOWN";
    }
}


MultiRootFindingMethod stringToMultiRootFindingMethod(const std::string& methodString) 
{
    // Iterate over the map with explicit type
    for (map<MultiRootFindingMethod, string>::const_iterator it = MultiRootFindingMethodMap.begin(); it != MultiRootFindingMethodMap.end(); ++it) 
    {
        if (it->second == methodString) 
        {
            return it->first;
        }
    }

    std::cout << "Invalid stringToMultiRootFindingMethod string: " << methodString << ". Aborting!" << std::endl;
    abort();
}


bool isValidMultiRootFindingMethod(const string& methodString)
{
    bool isMultiRootFindingMethodValid = false;
    // Iterate over the map with explicit type
    for (map<MultiRootFindingMethod, string>::const_iterator it = MultiRootFindingMethodMap.begin(); it != MultiRootFindingMethodMap.end(); ++it) 
    {
        if (it->second == methodString) 
        {
            isMultiRootFindingMethodValid = true;
            break;
        }
    }

    if( isMultiRootFindingMethodValid==false )
    {
        cout << "The value " + methodString + " is not a MultiRootFindingMethod!\n";
    }

    return isMultiRootFindingMethodValid;
}


//Multi-dimensional root-finding
void multiDimensionalRootFind(int n_eqs, double precision, double* x_init, void* params, int placeholder_f(const gsl_vector*, void*, gsl_vector*), MultiRootFindingMethod method)
{
	//set number of equations
	const size_t dim = n_eqs;

	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;
	
	//choose root fiding method
	if ( method==HYBRIDS )
	{
		T = gsl_multiroot_fsolver_hybrids;
	}
	else if ( method==HYBRID )
	{
		T = gsl_multiroot_fsolver_hybrid;
	}
	else if ( method==DNEWTON )
	{
		T = gsl_multiroot_fsolver_dnewton;
	}
	else if ( method==BROYDEN )
	{
		T = gsl_multiroot_fsolver_broyden;
	}
	else
    { 
        T = gsl_multiroot_fsolver_dnewton; //standart method if no other is selected
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

	while (status == GSL_CONTINUE && counter < MAX_ITERATIONS)
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


// One-dimensional root-finding
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
	else{ T = gsl_root_fsolver_bisection; }//standart method if no other is selected

	//allocate memory
	s =  gsl_root_fsolver_alloc (T);		

	gsl_root_fsolver_set (s, &F, x_low, x_high);

	//start iteration
	int counter = 0;
	int status = GSL_CONTINUE;

	double root;
	while (status == GSL_CONTINUE && counter < MAX_ITERATIONS){

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

	gsl_vector_free(xAux);
	gsl_vector_free(fAux);

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

    int status = gsl_multiroot_test_residual(fAux, precision);

    gsl_vector_free(xAux);
	gsl_vector_free(fAux);

    return status;
}


//sort a vector of complex numbers by absolute value
vector<gsl_complex> sortGSLComplexNumbersByAbsoluteSize(vector<gsl_complex> nonOrderedSet)
{	
	vector< tuple<double, int> > aux;
	for (int i = 0; i < int( nonOrderedSet.size() ); ++i)
	{
		double absoluteValue = gsl_complex_abs( nonOrderedSet[i] );
		aux.push_back( make_tuple(absoluteValue, i) );
	}
	sort(aux.begin(), aux.end());

	vector<gsl_complex> orderedSet;
	for (int i = 0; i < int( nonOrderedSet.size() ); ++i)
	{
		int index = get<1>(aux[i]);
		orderedSet.push_back( nonOrderedSet[index] );
	}

	if ( int(orderedSet.size())!=int(nonOrderedSet.size()) ){ cout << "Ordered set and non ordered set do not have the same size!\n"; abort(); }


	return orderedSet;
}


////////////////////////////////////////////////////////////////////////////////////////
//cubic equation solver using Cardano's method: a*x^3 + b*x^2 + c*x + d = 0

// term1_sqrt = -((b^2 c^2)/a^4), term2_sqrt = +((4 c^3)/a^3), term3_sqrt = +((4 b^3 d)/a^4), term4_sqrt = -((18 b c d)/a^3), term5_sqrt = +((27 d^2)/a^2)
// term1 = -((2 b^3)/a^3), term2 = +((9 b c)/a^2), term3 = -((27 d)/a)
gsl_complex cardanoA(gsl_complex a, gsl_complex b, gsl_complex c, gsl_complex d)
{   
    gsl_complex a2 = gsl_complex_pow_real(a, 2.0);   
    gsl_complex b2 = gsl_complex_pow_real(b, 2.0);
    gsl_complex c2 = gsl_complex_pow_real(c, 2.0);
    gsl_complex d2 = gsl_complex_pow_real(d, 2.0);
    gsl_complex a3 = gsl_complex_pow_real(a, 3.0);
    gsl_complex b3 = gsl_complex_pow_real(b, 3.0);
    gsl_complex c3 = gsl_complex_pow_real(c, 3.0);
    gsl_complex a4 = gsl_complex_pow_real(a, 4.0);
    gsl_complex bc = gsl_complex_mul(b, c);
    gsl_complex b2_c2 = gsl_complex_mul(b2, c2);
    gsl_complex b3d = gsl_complex_mul(b3, d);
    gsl_complex bcd = gsl_complex_mul(bc, d);

    //calculating terms inside the square root
    gsl_complex term1_sqrt = gsl_complex_mul_real(gsl_complex_div(b2_c2, a4), -1);
    gsl_complex term2_sqrt = gsl_complex_mul_real(gsl_complex_div(c3, a3), +4);
    gsl_complex term3_sqrt = gsl_complex_mul_real(gsl_complex_div(b3d, a4), +4);
    gsl_complex term4_sqrt = gsl_complex_mul_real(gsl_complex_div(bcd, a3), -18);
    gsl_complex term5_sqrt = gsl_complex_mul_real(gsl_complex_div(d2, a2), +27);

    //calculating the square root
    gsl_complex A_sqrt = gsl_complex_rect(0.0, 0.0);
    A_sqrt = gsl_complex_add(A_sqrt, term1_sqrt);
    A_sqrt = gsl_complex_add(A_sqrt, term2_sqrt);
    A_sqrt = gsl_complex_add(A_sqrt, term3_sqrt);
    A_sqrt = gsl_complex_add(A_sqrt, term4_sqrt);
    A_sqrt = gsl_complex_add(A_sqrt, term5_sqrt);
    A_sqrt = gsl_complex_sqrt(A_sqrt);

    gsl_complex term1 = gsl_complex_mul_real(gsl_complex_div(b3, a3), -2);
    gsl_complex term2 = gsl_complex_mul_real(gsl_complex_div(bc, a2), +9);
    gsl_complex term3 = gsl_complex_mul_real(gsl_complex_div(d, a), -27);

    gsl_complex A_aux = gsl_complex_rect(0.0, 0.0);
    A_aux = gsl_complex_add(A_aux, term1);
    A_aux = gsl_complex_add(A_aux, term2);
    A_aux = gsl_complex_add(A_aux, term3);
    A_aux = gsl_complex_add(A_aux, gsl_complex_mul_real(A_sqrt, 3.0*sqrt(3.0)) );
    A_aux = gsl_complex_pow_real(A_aux, 1./3);
    A_aux = gsl_complex_mul_real(A_aux, 1./(3.0*pow(2, 1./3)));

    return A_aux;
}


//B = (-9 a b^2 + 27 a^2 c)/(27 a^3)
gsl_complex cardanoB(gsl_complex a, gsl_complex b, gsl_complex c)
{   
    gsl_complex b2 = gsl_complex_pow_real(b, 2.0);
    gsl_complex ab2 = gsl_complex_mul(a, b2);
    gsl_complex a2 = gsl_complex_pow_real(a, 2.0);  
    gsl_complex a2c = gsl_complex_mul(a2, c);
    gsl_complex a3 = gsl_complex_pow_real(a, 3.0);

    gsl_complex B_aux = gsl_complex_rect(0.0, 0.0);
    B_aux = gsl_complex_mul_real(ab2, -9);
    B_aux = gsl_complex_add(B_aux, gsl_complex_mul_real(a2c, +27));
    B_aux = gsl_complex_mul_real(gsl_complex_div(B_aux, a3), 1./27);

    return B_aux;
}


vector<gsl_complex> solveCubicEquationCardano(gsl_complex a, gsl_complex b, gsl_complex c, gsl_complex d)
{   
    gsl_complex A = cardanoA(a, b, c, d);
    gsl_complex B = cardanoB(a, b, c);
    gsl_complex phi_m = gsl_complex_polar(1.0, -2.0*M_PI/3.0);
    gsl_complex phi_p = gsl_complex_polar(1.0, +2.0*M_PI/3.0);

    gsl_complex x1 = gsl_complex_rect(0.0, 0.0);
    x1 = gsl_complex_add(x1, gsl_complex_mul_real(gsl_complex_div(b, a), -1./3));
    x1 = gsl_complex_add(x1, A);
    x1 = gsl_complex_add(x1, gsl_complex_mul_real(gsl_complex_div(B, A), -1./3));

    gsl_complex x2 = gsl_complex_rect(0.0, 0.0);
    x2 = gsl_complex_add(x2, gsl_complex_mul_real(gsl_complex_div(b, a), -1./3));
    x2 = gsl_complex_add(x2, gsl_complex_mul(A,phi_m));
    x2 = gsl_complex_add(x2, gsl_complex_mul_real(gsl_complex_div(B, gsl_complex_mul(A,phi_m)), -1./3));

    gsl_complex x3 = gsl_complex_rect(0.0, 0.0);
    x3 = gsl_complex_add(x3, gsl_complex_mul_real(gsl_complex_div(b, a), -1./3));
    x3 = gsl_complex_add(x3, gsl_complex_mul(A,phi_p));
    x3 = gsl_complex_add(x3, gsl_complex_mul_real(gsl_complex_div(B, gsl_complex_mul(A,phi_p)), -1./3));

    vector<gsl_complex> solutions;
    solutions.push_back( x1 );
    solutions.push_back( x2 );
    solutions.push_back( x3 );
    solutions = sortGSLComplexNumbersByAbsoluteSize(solutions);

    return solutions;
}


vector<gsl_complex> calculateEigenvalues3By3ComplexMatrix(ComplexSquareMatrixGSL M)
{	
	if ( int(M.getDimension())!=3 ){ cout << "Trying to calculate eigenvalues of matrix with size !=3 using Cardano!\n"; abort(); }

	//write cubic polynomial for eigenvalues
	gsl_complex a = gsl_complex_rect(-1.0, 0.0);

    gsl_complex b = gsl_complex_rect(0.0, 0.0);
    b = gsl_complex_add(b, M.getValue(0, 0));
    b = gsl_complex_add(b, M.getValue(1, 1));
    b = gsl_complex_add(b, M.getValue(2, 2));
    
    gsl_complex c = gsl_complex_rect(0.0, 0.0);
    c = gsl_complex_add(c, gsl_complex_mul( M.getValue(0, 1) , M.getValue(1, 0) ));
    c = gsl_complex_sub(c, gsl_complex_mul( M.getValue(0, 0) , M.getValue(1, 1) ));
    c = gsl_complex_add(c, gsl_complex_mul( M.getValue(0, 2) , M.getValue(2, 0) ));
    c = gsl_complex_add(c, gsl_complex_mul( M.getValue(1, 2) , M.getValue(2, 1) ));
    c = gsl_complex_sub(c, gsl_complex_mul( M.getValue(0, 0) , M.getValue(2, 2) ));
    c = gsl_complex_sub(c, gsl_complex_mul( M.getValue(1, 1) , M.getValue(2, 2) ));

    gsl_complex d = gsl_complex_rect(0.0, 0.0);
    d = gsl_complex_sub(d, gsl_complex_mul( M.getValue(0, 2) , gsl_complex_mul( M.getValue(2, 0) , M.getValue(1, 1) ) ));
    d = gsl_complex_add(d, gsl_complex_mul( M.getValue(2, 0) , gsl_complex_mul( M.getValue(0, 1) , M.getValue(1, 2) ) ));
    d = gsl_complex_add(d, gsl_complex_mul( M.getValue(0, 2) , gsl_complex_mul( M.getValue(2, 1) , M.getValue(1, 0) ) ));
    d = gsl_complex_sub(d, gsl_complex_mul( M.getValue(0, 0) , gsl_complex_mul( M.getValue(1, 2) , M.getValue(2, 1) ) ));
    d = gsl_complex_sub(d, gsl_complex_mul( M.getValue(0, 1) , gsl_complex_mul( M.getValue(1, 0) , M.getValue(2, 2) ) ));
    d = gsl_complex_add(d, gsl_complex_mul( M.getValue(0, 0) , gsl_complex_mul( M.getValue(1, 1) , M.getValue(2, 2) ) ));

    //use cardano to solve the cubic equation
    vector<gsl_complex> eigenvalues = solveCubicEquationCardano(a, b, c, d);

    return eigenvalues;
}


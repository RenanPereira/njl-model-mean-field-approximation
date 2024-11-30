#include "gsl_wrapper/Interpolation1DimGSL.h"


void Interpolation1DimGSL::setSpline()
{   
	const int N = interpolationSize;

    //create arrays necessary for interpolation
    double variable[N];
    copy(discretizedVariable.begin(), discretizedVariable.end(), variable);
        
    double function[N];
    copy(discretizedFunction.begin(), discretizedFunction.end(), function);

    //allocate memory for GSL interpolation
    if ( method==linear ){ spline = gsl_spline_alloc(gsl_interp_linear, N); }
    else if( method==steffen ){ spline = gsl_spline_alloc(gsl_interp_steffen, N); }
    else if( method==cubic ){ spline = gsl_spline_alloc(gsl_interp_cspline, N); }
    else if( method==akima ){ spline = gsl_spline_alloc(gsl_interp_akima, N); }	
    else if( method==polynomial ){ spline = gsl_spline_alloc(gsl_interp_polynomial, N); }
    else { spline = gsl_spline_alloc(gsl_interp_linear, N); }

    //initiate interpolation
    gsl_spline_init(spline, variable, function, N);
}


void Interpolation1DimGSL::unsetSpline()
{   
	gsl_spline_free(spline);//clean memory
}


void Interpolation1DimGSL::setAccelerator()
{   
	accelerator = gsl_interp_accel_alloc();
}


void Interpolation1DimGSL::unsetAccelerator()
{   
	gsl_interp_accel_free(accelerator);//clean memory
}


Interpolation1DimGSL::Interpolation1DimGSL(InterpolationMethod methodAux, vector<double> discretizedVariableAux, vector<double> discretizedFunctionAux)
{   
    method = methodAux;
    discretizedVariable = discretizedVariableAux;
    discretizedFunction = discretizedFunctionAux;
    interpolationSize = discretizedFunction.size();
    interpolationLowerBound = discretizedVariable[0];
    interpolationUpperBound = discretizedVariable[interpolationSize-1];
        
    //check if the provided lists discretizedVariable and discretizedFunction have the same size
    if ( interpolationSize!=int(discretizedVariable.size()) )
    {
        cout << "The list x and y values, provided for interpolation, does not have the same size!" << "\n"; 
        abort();
    }

    //check if the independent variable is always increasing: the GSL libraray requires it!
    for (int i = 0; i < interpolationSize-1; i++)
    {
        if( discretizedVariable[i+1]<discretizedVariable[i] )
        {
            cout << "The independent variable is not always increasing, GSL will not know what to do!\n";
        }
    }

    setSpline();
    setAccelerator();
}


Interpolation1DimGSL::Interpolation1DimGSL(InterpolationMethod methodAux, double xMin, double xMax, int N, OneVariableFunction function)
{   
	vector<double> discretizedVariableAux(N,0.0);
	vector<double> discretizedFunctionAux(N,0.0);

	double delta = (xMax - xMin)/( N - 1 );
	for (int i = 0; i < N; i++)
	{	
		double x = xMin + i*delta;
		discretizedVariableAux[i] = x;
		discretizedFunctionAux[i] = function.evaluate(x);
	}

	method = methodAux;
    discretizedVariable = discretizedVariableAux;
    discretizedFunction = discretizedFunctionAux;
    interpolationSize = discretizedFunction.size();
    interpolationLowerBound = discretizedVariable[0];
    interpolationUpperBound = discretizedVariable[interpolationSize-1];
        
    //check if the provided lists discretizedVariable and discretizedFunction have the same size
    if ( interpolationSize!=int(discretizedVariable.size()) )
    {
        cout << "The list x and y values, provided for interpolation, does not have the same size!" << "\n"; 
        abort();
    }

    //check if the independent variable is always increasing: the GSL libraray requires it!
    for (int i = 0; i < interpolationSize-1; i++)
    {
        if( discretizedVariable[i+1]<discretizedVariable[i] )
        {
            cout << "The independent variable is not always increasing, GSL will not know what to do!\n";
        }
    }

    setSpline();
    setAccelerator();
}


bool Interpolation1DimGSL::valueInInterpolationInterval(double x)
{   
    double xMin = getInterpolationLowerBound();
	double xMax = getInterpolationUpperBound();
	bool test = (x>=xMin && x<=xMax);

	return test;
}


double Interpolation1DimGSL::evaluate(double x)
{	
	double fx = 0.0;
	if ( valueInInterpolationInterval(x) )
	{	
		//find the value of the interpolated function at some specific x-value
		fx = gsl_spline_eval(spline, x, accelerator);
	}
	else
	{ 
		cout << "Evaluating the interpolation outside of the interpolation range!\n";
		fx = 1.0/0.0; 
	}

	return fx;
}


double Interpolation1DimGSL::evaluate1stDerivative(double x)
{
	double fx = 0.0;
	if ( valueInInterpolationInterval(x) )
	{	
		//find the value of the interpolated function at some specific x-value
		fx = gsl_spline_eval_deriv(spline, x, accelerator);
	}
	else
	{ 	
		cout << "Evaluating 1st derivative of the interpolation outside of the interpolation range!\n";
		fx = 1.0/0.0; 
	}

	return fx;
}


double Interpolation1DimGSL::evaluate2ndDerivative(double x)
{
	double fx = 0.0;
	if ( valueInInterpolationInterval(x) )
	{	
		//find the value of the interpolated function at some specific x-value
		fx = gsl_spline_eval_deriv2(spline, x, accelerator);
	}
	else
	{ 	
		cout << "Evaluating 2nd derivative of the interpolation outside of the interpolation range!\n";
		fx = 1.0/0.0; 
	}

	return fx;
}


double Interpolation1DimGSL::evaluateIntegral(double a, double b)
{	
	double integral = 0.0;
	if ( valueInInterpolationInterval(a) && valueInInterpolationInterval(b) )
	{	
		integral = gsl_spline_eval_integ(spline, a, b, accelerator);
	}
	else
	{ 	
		cout << "The integration interval is outside of the interpolation range!\n";
		integral = 1.0/0.0; 
	}

	return integral;
}


vector<double> Interpolation1DimGSL::findRoots(RootFindingMethod methodAux, double precision)
{	
	const int N = interpolationSize;
	
	vector<int> guesses;
	for (int i = 0; i < N-1; ++i)
	{
		if      ( discretizedFunction[i] < 0 && discretizedFunction[i+1]>0 ){ guesses.push_back(i); }
		else if ( discretizedFunction[i] > 0 && discretizedFunction[i+1]<0 ){ guesses.push_back(i); }
	}

	vector<double> roots;
	if ( guesses.size()!=0 )
	{
		for (int i = 0; i < int(guesses.size()); ++i)
		{	
			double prevGuess = discretizedVariable[guesses[i]];
			double nextGuess = discretizedVariable[guesses[i]+1];
			double roots_aux = OneDimensionalRootFind(precision, prevGuess, nextGuess, this, &equationToFindRootsOfInterpolation, methodAux);
			roots.push_back( roots_aux );
		}
	}

	return roots;
}


vector<double> Interpolation1DimGSL::findRoots1stDerivative(RootFindingMethod methodAux, double precision)
{	
	const int N = interpolationSize;
	
	vector<int> guesses;
	for (int i = 0; i < N-1; ++i)
	{	
		double xi = discretizedVariable[i];
		double xiPlusOne = discretizedVariable[i+1];
		double discretizedFunction1stDerivativeAtXi = evaluate1stDerivative(xi);
		double discretizedFunction1stDerivativeAtPlusOne = evaluate1stDerivative(xiPlusOne);
		if      ( discretizedFunction1stDerivativeAtXi < 0 && discretizedFunction1stDerivativeAtPlusOne>0 ){ guesses.push_back(i); }
		else if ( discretizedFunction1stDerivativeAtXi > 0 && discretizedFunction1stDerivativeAtPlusOne<0 ){ guesses.push_back(i); }
	}

	vector<double> roots;
	if ( guesses.size()!=0 )
	{
		for (int i = 0; i < int(guesses.size()); ++i)
		{	
			double prevGuess = discretizedVariable[guesses[i]];
			double nextGuess = discretizedVariable[guesses[i]+1];
			double roots_aux = OneDimensionalRootFind(precision, prevGuess, nextGuess, this, &equationToFindRootsOfInterpolation1stDerivative, methodAux);
			roots.push_back( roots_aux );
		}
	}

	return roots;
}


vector<double> Interpolation1DimGSL::findRoots2ndDerivative(RootFindingMethod methodAux, double precision)
{	
	const int N = interpolationSize;
	
	vector<int> guesses;
	for (int i = 1; i < N-1; ++i)
	{	
		double xi = discretizedVariable[i];
		double xiPlusOne = discretizedVariable[i+1];
		double discretizedFunction2ndDerivativeAtXi = evaluate2ndDerivative(xi);
		double discretizedFunction2ndDerivativeAtPlusOne = evaluate2ndDerivative(xiPlusOne);
		if      ( discretizedFunction2ndDerivativeAtXi < 0 && discretizedFunction2ndDerivativeAtPlusOne>0 ){ guesses.push_back(i); }
		else if ( discretizedFunction2ndDerivativeAtXi > 0 && discretizedFunction2ndDerivativeAtPlusOne<0 ){ guesses.push_back(i); }
	}

	vector<double> roots;
	if ( guesses.size()!=0 )
	{
		for (int i = 0; i < int(guesses.size()); ++i)
		{	
			double prevGuess = discretizedVariable[guesses[i]];
			double nextGuess = discretizedVariable[guesses[i]+1];
			double roots_aux = OneDimensionalRootFind(precision, prevGuess, nextGuess, this, &equationToFindRootsOfInterpolation2ndDerivative, methodAux);
			roots.push_back( roots_aux );
		}
	}

	return roots;
}


void Interpolation1DimGSL::tests(OneVariableFunction testFunction)
{
	vector<double> roots = findRoots(brent, 1E-8);

	cout << "Number of simple roots found = " << roots.size() << "\n";
	cout << "The roots in the provided interpolation bounds are (x , y):\n";
	for (int i = 0; i < int(roots.size()); i++)
	{	
		double x = roots[i];
		cout << "(" << x << " , " << testFunction.evaluate(x) << ")" << "\n";
	}

    vector<double> extrema = findRoots1stDerivative(brent, 1E-8);
    cout << "Number of extrema found = " << extrema.size() << "\n";
	cout << "The extrema in the provided interpolation bounds are  (x , y):\n";
    for (int i = 0; i < int(extrema.size()); i++)
    {	
		double x = extrema[i];
        cout << "(" << x << " , " << testFunction.evaluate(x) << ")" << "\n";
    }

	double A = getInterpolationLowerBound();
	double B = getInterpolationUpperBound();
	double integral = evaluateIntegral(A, B);
	cout << "Integral of the interpolation in the entire interpolation range: " << integral << "\n";
}


double polynomialTest(double x, void* parameters)
{	
	(void)(parameters); /* avoid unused parameter warning */
	//double poly = pow(x,6) - 5*pow(x,5) - 12*pow(x,4) - 3*pow(x,3) - 24*pow(x,2) + 7*pow(x,1) + 6;
	double poly = (x-1)*(x-2)*(x-3)*(x-4)*(x+1)*(x+2)*(x+3)*(x+4)/1000.0;

	return poly;
}


double equationToFindRootsOfInterpolation(double x, void *interpolationAux)
{
	//InterpolationOfDiscreteFunction1D *interpolation = static_cast<InterpolationOfDiscreteFunction1D *>(interpolationAux);
	Interpolation1DimGSL* interpolation = ((class Interpolation1DimGSL *)(interpolationAux));

	double fx = interpolation->evaluate(x);

	return fx;
}


double equationToFindRootsOfInterpolation1stDerivative(double x, void *interpolationAux)
{
	Interpolation1DimGSL* interpolation = ((class Interpolation1DimGSL *)(interpolationAux));

	double fx = interpolation->evaluate1stDerivative(x);

	return fx;
}


double equationToFindRootsOfInterpolation2ndDerivative(double x, void *interpolationAux)
{
	Interpolation1DimGSL* interpolation = ((class Interpolation1DimGSL *)(interpolationAux));

	double fx = interpolation->evaluate2ndDerivative(x);

	return fx;
}




#ifndef INTEGRATION1DIMGSL_H
#define INTEGRATION1DIMGSL_H


#include <gsl/gsl_integration.h>
#include <vector>

using namespace std;


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


class Integration1DimGSL
{	
public:
	double lowerBound;
	double upperBound;
	double minimumDistanceBetweenBounds = 1E-14;

	GeneralIntegrandParameters* integrandParameters;
	gsl_function F;

	double absolutePrecision;
	double relativePrecision;

	size_t workspaceLimitSize;

	double result;
	double error;
	size_t numberFunctionEvaluations;

public:
	void setVariables(double , double , GeneralIntegrandParameters* , double integrand(double, void*), double , double , int );

	double getLowerBound(){ return lowerBound; }
	double getUpperBound(){ return upperBound; }

	double getAbsolutePrecision(){ return absolutePrecision; }
	double getRelativePrecision(){ return relativePrecision; }

	GeneralIntegrandParameters* getIntegrandParameters(){ return integrandParameters; }

	double getResult(){ return result; }
	double getError(){ return error; }
	int getNumberEvaluations(){ return int(numberFunctionEvaluations); }

	virtual double evaluate();

	void errorHandler(int , string );
};



/*
The QNG algorithm is a non-adaptive procedure which uses fixed Gauss-Kronrod-Patterson abscissae to sample the
integrand at a maximum of 87 points. It is provided for fast integration of smooth functions.
This function applies the Gauss-Kronrod 10-point, 21-point, 43-point and 87-point integration rules in succession.
*/
class Integration1DimGSLQNG : public Integration1DimGSL
{	
public:
	Integration1DimGSLQNG(double , double , GeneralIntegrandParameters* , double (double, void*), double , double );

	double evaluate() override;
};



/*
The QAG algorithm is a simple adaptive integration procedure. The integration region is divided into subintervals, and
on each iteration the subinterval with the largest estimated error is bisected. This reduces the overall error rapidly, as
the subintervals become concentrated around local difficulties in the integrand. 
The integration rule is determined by the value of key: 1-6, corresponding to the 15, 21, 31, 41, 51 and 61 point Gauss-Kronrod rules. 
The higher-order rules give better accuracy for smooth functions, while lower-order rules save time when the function contains 
local difficulties, such as discontinuities.
On each iteration the adaptive integration strategy bisects the interval with the largest error estimate. The subin-
tervals and their results are stored in the memory provided by workspace. The maximum number of subintervals
is given by limit, which may not exceed the allocated size of the workspace.
*/
class Integration1DimGSLQAG : public Integration1DimGSL
{	
public:
	int key;

public:
	Integration1DimGSLQAG(double , double , GeneralIntegrandParameters* , double (double, void*), double , double , int, int );

	double evaluate() override;
};


/*
The presence of an integrable singularity in the integration region causes an adaptive routine to concentrate new 
subintervals around the singularity. As the subintervals decrease in size the successive approximations to the integral 
converge  in a limiting fashion. This approach to the limit can be accelerated using an extrapolation procedure. The QAGS
algorithm combines adaptive bisection with the Wynn epsilon-algorithm to speed up the integration of many types of
integrable singularities.
This function applies the Gauss-Kronrod 21-point integration rule adaptively. The results are extrapolated using the 
epsilon-algorithm, which accelerates the convergence of the integral in the presence of
discontinuities and integrable singularities.
*/
class Integration1DimGSLQAGS : public Integration1DimGSL
{	
public:
	Integration1DimGSLQAGS(double , double , GeneralIntegrandParameters* , double (double, void*), double , double , int );

	double evaluate() override;
};


/*
This function computes the Cauchy principal value of the integral.
The adaptive bisection algorithm of QAG is used, with modifications to ensure that subdivisions do not occur at
the singular point. When a subinterval contains the singular point or is close to it then a special 25-point
modified Clenshaw-Curtis rule is used to control the singularity. Further away from the singularity the algorithm
uses an ordinary 15-point Gauss-Kronrod integration rule.
*/
class Integration1DimGSLQAWC : public Integration1DimGSL
{	
public:
	double singularity;

public:
	Integration1DimGSLQAWC(double , double , double, GeneralIntegrandParameters* , double (double, void*), double , double , int );

	double evaluate() override;
};


/*
This function applies the adaptive integration algorithm QAGS taking account of the user-supplied locations of
singular points. The vector singularities should contain the endpoints of the integration ranges defined
by the integration region and locations of the singularities. 
*/
class Integration1DimGSLQAGP : public Integration1DimGSL
{	
public:
	vector<double> singularities;

public:
	Integration1DimGSLQAGP(double , double , vector<double>, GeneralIntegrandParameters* , double (double, void*), double , double , int );

	double evaluate() override;
};


/*
CQUAD is a new doubly-adaptive general-purpose quadrature routine which can handle most types of singularities,
non-numerical function values such as Inf or NaN, as well as some divergent integrals.
It generally requires more function evaluations than the integration routines in QUADPACK, yet fails less often for difficult integrands.
The underlying algorithm uses a doubly-adaptive scheme in which Clenshaw-Curtis quadrature rules of increasing
degree are used to compute the integral in each interval. The L2-norm of the difference between the underlying inter-
polatory polynomials of two successive rules is used as an error estimate. The interval is subdivided if the difference
between two successive rules is too large or a rule of maximum degree has been reached.
*/
class Integration1DimGSLCQUAD : public Integration1DimGSL
{	
public:
	Integration1DimGSLCQUAD(double , double , GeneralIntegrandParameters* , double (double, void*), double , double , int );

	double evaluate() override;
};


/*
This function computes the integral of the function f over the infinite interval (−∞, +∞).
The integral is mapped onto the semi-open interval (0, 1] using the transformation x = (1 − t)/t.
It is then integrated using the QAGS algorithm. The normal 21-point Gauss-Kronrod rule of QAGS is replaced
by a 15-point rule, because the transformation can generate an integrable singularity at the origin. In this case a
lower-order rule is more efficient.
*/
class Integration1DimGSLQAGI : public Integration1DimGSL
{	
public:
	Integration1DimGSLQAGI(GeneralIntegrandParameters* , double (double, void*), double , double , int );

	double evaluate() override;
};


/*
This function computes the integral of the function f over the semi-infinite interval (a, +∞). The integral is
mapped onto the semi-open interval (0, 1] using the transformation x = a + (1 − t)/t,
and then integrated using the QAGS algorithm.
*/
class Integration1DimGSLQAGIU : public Integration1DimGSL
{	
public:
	Integration1DimGSLQAGIU(double , GeneralIntegrandParameters* , double (double, void*), double , double , int );

	double evaluate() override;
};


/*
This function computes the integral of the function f over the semi-infinite interval (−∞, b). The integral is
mapped onto the semi-open interval (0, 1] using the transformation x = b − (1 − t)/t,
and then integrated using the QAGS algorithm.
*/
class Integration1DimGSLQAGIL : public Integration1DimGSL
{	
public:
	Integration1DimGSLQAGIL(double , GeneralIntegrandParameters* , double (double, void*), double , double , int );

	double evaluate() override;
};


/*
The QAWS algorithm is designed for integrands with algebraic-logarithmic singularities at the end-points of an integration region. 
In order to work efficiently the algorithm requires a precomputed table of Chebyshev moments.
This function allocates space for a gsl_integration_qaws_tablestruct describing a singular weight function
w(x) with the parameters (α, β, µ, ν),
w(x) = (x − a)α(b − x)β logµ(x − a) logν(b − x)
where α > −1, β > −1, and µ = 0, 1, ν = 0, 1. The weight function can take four different forms depending
on the values of µ and ν.
The adaptive bisection algorithm of QAG is used. When a subinterval contains one of the endpoints then a
special 25-point modified Clenshaw-Curtis rule is used to control the singularities. For subintervals which do
not include the endpoints an ordinary 15-point Gauss-Kronrod integration rule is used.
*/
class Integration1DimGSLQAWS : public Integration1DimGSL
{	
public:
	double alpha;
	double beta;
	int mu;
	int nu;

public:
	Integration1DimGSLQAWS(double , double , GeneralIntegrandParameters* , double (double, void*), double , double , int , double , double , int , int );

	double evaluate() override;
};


double changeQAWCIntegrandToQAGSIntegrand(double , void *);

class ChangeQAWCToQAGSParameters : public GeneralIntegrandParameters
{	
private:
	GeneralIntegrandParameters* numeratorParameters;
	double (* numeratorQAWC)(double, void * parameters);
	double singularity;

public:
	ChangeQAWCToQAGSParameters(){};
	ChangeQAWCToQAGSParameters(GeneralIntegrandParameters* numeratorParametersAux, double numeratorQAWCAux(double, void*), double singularityAux)
	{	
	    numeratorParameters = numeratorParametersAux;
	    numeratorQAWC = numeratorQAWCAux;
    	singularity = singularityAux;
	};

	ChangeQAWCToQAGSParameters(void* auxiliar)
    {   
        numeratorParameters = ((class ChangeQAWCToQAGSParameters *)(auxiliar))->numeratorParameters;
        numeratorQAWC = ((class ChangeQAWCToQAGSParameters *)(auxiliar))->numeratorQAWC;
        singularity = ((class ChangeQAWCToQAGSParameters *)(auxiliar))->singularity;
    };

    double getSingularity(){ return singularity; };
    double evaluateNumerator(double x){ double numerator = numeratorQAWC(x, numeratorParameters); return numerator; };

	void printIntegrandVariables() override
    {   
    	numeratorParameters->printIntegrandVariables();
        cout << "singularity = " << singularity << "\n";
    }
};


class Integration1DimGSLQAWCQAGS : public Integration1DimGSL
{
private:
	ChangeQAWCToQAGSParameters numeratorParameters;
	double singularity;
	double minimumDistanceBetweenBoundAndSingularity = 0.0;

public:
	Integration1DimGSLQAWCQAGS(double , double , double, GeneralIntegrandParameters* , double (double, void*), double , double , int );
	Integration1DimGSLQAWCQAGS(double , double , double, GeneralIntegrandParameters* , double (double, void*), double , double , int , double);

	bool isSingularityInsideTheIntegrationInterval();

	double getSingularity(){ return singularity; }

	double evaluate() override;
};


////////////////////////////////////////////////////////////////////////////////////////


enum TrapezoidalRule { normal, alternative };

class CompositeTrapezoidalSum
{	
public:
	double lowerBound;
	double upperBound;
	int numberOfPartitions;

	GeneralIntegrandParameters* integrandParameters;
	double (* integrand)(double, void * parameters);

	TrapezoidalRule rule;
	
	double result;

public:
	CompositeTrapezoidalSum(double , double , int , GeneralIntegrandParameters* , double (double, void*) );
	CompositeTrapezoidalSum(double , double , int , GeneralIntegrandParameters* , double (double, void*) , TrapezoidalRule );

	void setVariables(double , double , int , GeneralIntegrandParameters* , double (double, void*), TrapezoidalRule );

	double evaluateNormal();
	double evaluateAlternative();
	double evaluate();
	double evaluateCauchyPV(double singularity);
};


////////////////////////////////////////////////////////////////////////////////////////


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

double integrandTest(double , void *);

double integrandTestCauchy(double , void *);

double integrandTestQAGP(double , void *);

double integrandTestQAGI(double , void *);

double integrandRiemannCPV(double , void *);

double integrandTestQAWS(double , void *);

void testIntegration1DimGSL();

void testCompositeTrapezoidalSum();


#endif

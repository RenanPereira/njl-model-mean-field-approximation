#include <iostream>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include "Integration1DimGSL.h"

using namespace std;


void Integration1DimNewtonCotes::setVariables(double lowerBoundAux, double upperBoundAux, int numberOfPartitionsAux, GeneralIntegrandParameters* integrandParametersAux, double integrandAux(double, void*), NewtonCotesRule ruleAux)
{
    lowerBound = lowerBoundAux;
    upperBound = upperBoundAux;
    numberOfPartitions = numberOfPartitionsAux;

    integrandParameters = integrandParametersAux;
    integrand = integrandAux;

    rule = ruleAux;

    //check if number of partitions are appropriate for the chosen rule
    if ( rule==trapezoidal )
    {
        if ( numberOfPartitions<2 )
        {
            cout << "Integration1DimNewtonCotes: to use the trapezoidal rule, at least 2 partitions are necessary!\n";
            abort();
        }
    }

    if ( rule==alternativeCompositeSimpson )
    {
        if ( numberOfPartitions<8 )
        {
            cout << "Integration1DimNewtonCotes: to use the alternative Composite Simpson rule, at least 8 partitions are necessary!\n";
            abort();
        }
    }
}


Integration1DimNewtonCotes::Integration1DimNewtonCotes(double lowerBoundAux, double upperBoundAux, int numberOfPartitionsAux, GeneralIntegrandParameters* integrandParametersAux, double integrandAux(double, void*))
{
    setVariables(lowerBoundAux, upperBoundAux, numberOfPartitionsAux, integrandParametersAux, integrandAux, trapezoidal);
}


Integration1DimNewtonCotes::Integration1DimNewtonCotes(double lowerBoundAux, double upperBoundAux, int numberOfPartitionsAux, GeneralIntegrandParameters* integrandParametersAux, double integrandAux(double, void*), NewtonCotesRule ruleAux)
{
    setVariables(lowerBoundAux, upperBoundAux, numberOfPartitionsAux, integrandParametersAux, integrandAux, ruleAux);
}


double Integration1DimNewtonCotes::evaluateTrapezoidal()
{   
    
    double dx = (upperBound-lowerBound)/(numberOfPartitions-1);
    double area = 0.0;
    for (int i = 1; i < numberOfPartitions-1; ++i)
    {   
        area = area + integrand(lowerBound + i*dx, integrandParameters);
    }

    area = dx*( area + 0.5*integrand(lowerBound, integrandParameters) + 0.5*integrand(upperBound, integrandParameters) );
    result = area;
    
    return area;
}


double Integration1DimNewtonCotes::evaluateAlternativeCompositeSimpson()
{   
    double dx = (upperBound-lowerBound)/(numberOfPartitions-1);
    double area = 0.0;
    for (int i = 4; i < numberOfPartitions-4; ++i)
    {   
        area = area + integrand(lowerBound + i*dx, integrandParameters);
    }

    area = ( dx/(48.0) )*( + 17.0*integrand(lowerBound + 0.0*dx, integrandParameters)
                           + 59.0*integrand(lowerBound + 1.0*dx, integrandParameters)
                           + 43.0*integrand(lowerBound + 2.0*dx, integrandParameters)
                           + 49.0*integrand(lowerBound + 3.0*dx, integrandParameters)
                           + 48.0*area
                           + 49.0*integrand(upperBound - 3.0*dx, integrandParameters)
                           + 43.0*integrand(upperBound - 2.0*dx, integrandParameters)
                           + 59.0*integrand(upperBound - 1.0*dx, integrandParameters)
                           + 17.0*integrand(upperBound - 0.0*dx, integrandParameters)
                         );
    result = area;

    return area;
}


double Integration1DimNewtonCotes::evaluate()
{   
    double area = 0.0;

    if ( rule==trapezoidal )
    {
        area = evaluateTrapezoidal();
    }
    else if ( rule==alternativeCompositeSimpson )
    {
        area = evaluateAlternativeCompositeSimpson();
    }

    return area;
}


double Integration1DimNewtonCotes::evaluateAvoidingSingularPoint(double singularity)
{ 

    double dx = (upperBound-lowerBound)/(numberOfPartitions-1);
    double delta = 2*dx;

    double area = 0.0;

    double a = 0;
    double b = 0;
    double orientation = +1;
    if ( lowerBound<upperBound )
    {
        a = lowerBound;
        b = upperBound;
        orientation = +1;
    }
    else
    {
        a = upperBound;
        b = lowerBound;
        orientation = -1;   
    }

    if ( fabs(b-singularity)<fabs(a-singularity) )
    {
        Integration1DimNewtonCotes LowRiemannSum(a, singularity-delta - fabs(b - (singularity+delta)), numberOfPartitions, integrandParameters, integrand, rule);
        area = area + LowRiemannSum.evaluate();
    
        Integration1DimNewtonCotes MiddleRiemannSum2(singularity-delta - fabs(b - (singularity+delta)), singularity-delta, numberOfPartitions, integrandParameters, integrand, rule);
        area = area + MiddleRiemannSum2.evaluate();

        Integration1DimNewtonCotes UpperRiemannSum(singularity+delta, b, numberOfPartitions, integrandParameters, integrand, rule);
        area = area + UpperRiemannSum.evaluate();

        area = orientation*area;
    }
    else
    {
        Integration1DimNewtonCotes LowRiemannSum(a, singularity-delta, numberOfPartitions, integrandParameters, integrand, rule);
        area = area + LowRiemannSum.evaluate();
        
        Integration1DimNewtonCotes MiddleRiemannSum2(singularity+delta, (singularity+delta) + fabs(a - (singularity+delta)), numberOfPartitions, integrandParameters, integrand, rule);
        area = area + MiddleRiemannSum2.evaluate();

        Integration1DimNewtonCotes UpperRiemannSum((singularity+delta) + fabs(a - (singularity+delta)), b, numberOfPartitions, integrandParameters, integrand, rule);
        area = area + UpperRiemannSum.evaluate();
    }
    
    return area;

/*
    double dx = (upperBound-lowerBound)/(numberOfPartitions-1);
    double delta = 2*dx;

    double area = 0.0;
    
    Integration1DimNewtonCotes LowRiemannSum(lowerBound, singularity-delta, numberOfPartitions, integrandParameters, integrand, rule);
    area = area + LowRiemannSum.evaluate();
    
    Integration1DimNewtonCotes UpperRiemannSum(singularity+delta, upperBound, numberOfPartitions, integrandParameters, integrand, rule);
    area = area + UpperRiemannSum.evaluate();
    
    return area;
*/
}


////////////////////////////////////////////////////////////////////////////////////////


void Integration1DimGSL::setVariables(double lowerBoundAux, double upperBoundAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux)
{
    lowerBound = lowerBoundAux;
    upperBound = upperBoundAux;
    integrandParameters = integrandParametersAux;

    F.function = integrand;
    F.params = integrandParameters;

    absolutePrecision = absolutePrecisionAux;
    relativePrecision = relativePrecisionAux;


    //check workspace
    if ( workspaceLimitSizeAux<3 ){ printf("Problem with chosen workspace size! It must be larger than 2! \n"); abort(); }
    workspaceLimitSize = workspaceLimitSizeAux;   
}


double Integration1DimGSL::evaluate()
{
    cout << "Using evaluate from Integration1DimGSL! Such method is not defined!\n";

    return 0;
}


void Integration1DimGSL::errorHandler(int code, string methodName)
{
    if( code!=0 )
    {   
        //print the integral method that failed
        cout << "\n";
        cout << "Problem in the integration using the " + methodName + " method! " << "\n";

        //print the type of error
        if      (code==GSL_EMAXITER){ cout << gsl_strerror(GSL_EMAXITER) << "\n"; }
        else if (code==GSL_EDIVERGE){ cout << gsl_strerror(GSL_EDIVERGE) << "\n"; }
        else if (code==GSL_ESING)   { cout << gsl_strerror(GSL_ESING)    << "\n"; }
        else if (code==GSL_EROUND)  { cout << gsl_strerror(GSL_EROUND)   << "\n"; }
        else if (code==GSL_EDOM)    { cout << gsl_strerror(GSL_EDOM)     << "\n"; }
        else                        { cout << "Error code not identified in the GSL library!\n"; }

        cout << "\n";
        cout << "Integration1DimGSL variables:\n";
        cout << "lower bound = " << lowerBound << "\n";
        cout << "upper bound = " << upperBound << "\n";
        cout << "absolute precision = "  << absolutePrecision << "\n";
        cout << "relative precision = "  << relativePrecision << "\n";
        cout << "\n";

        //print the integrand parameters used in the failed integration for debugging
        integrandParameters->printIntegrandVariables();

        //call the method for the behaviour after failed integration
        integrandParameters->behaviourAfterFailedIntegration();
    }
}


//QNG method constructor
Integration1DimGSLQNG::Integration1DimGSLQNG(double lowerBoundAux, double upperBoundAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux)
{   
    setVariables(lowerBoundAux, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, 3);
}


//QNG method integral evaluator
double Integration1DimGSLQNG::evaluate()
{   
    if ( fabs(lowerBound-upperBound)>minimumDistanceBetweenBounds )
    {
        //save original handler, turn off the error handler
        gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


        int code = gsl_integration_qng(&F, lowerBound, upperBound, absolutePrecision, relativePrecision, &result, &error, &numberFunctionEvaluations);
        
        //handle possible integration error
        errorHandler(code, "QNG");


        //reset the error handler to previous state before switching off
        gsl_set_error_handler(old_error_handler);
    }
    else
    {
        result = 0.0;
    }

    return result;
}


//QAG method constructor
Integration1DimGSLQAG::Integration1DimGSLQAG(double lowerBoundAux, double upperBoundAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux, int keyAux)
{   
    setVariables(lowerBoundAux, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);

    //check key
    if ( keyAux<1 || keyAux>6 ){ printf("Using QAG method: problem with chosen key! Value must be an integer in the interval [1,6]! \n"); abort(); }
    key = keyAux;
}


//QAG method integral evaluator
double Integration1DimGSLQAG::evaluate()
{   
    if ( fabs(lowerBound-upperBound)>minimumDistanceBetweenBounds )
    {
        //save original handler, turn off the error handler
        gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


        //allocate memory
        gsl_integration_workspace *workPtr = gsl_integration_workspace_alloc(workspaceLimitSize);
        
        int code = gsl_integration_qag(&F, lowerBound, upperBound, absolutePrecision, relativePrecision, workspaceLimitSize, key, workPtr, &result, &error);

        //handle possible integration error
        errorHandler(code, "QAG");

        //clean memory
        gsl_integration_workspace_free(workPtr);

        
        //reset the error handler to previous state before switching off
        gsl_set_error_handler(old_error_handler);
    }
    else
    {
        result = 0.0;
    }

    return result;
}


//QAGS method constructor
Integration1DimGSLQAGS::Integration1DimGSLQAGS(double lowerBoundAux, double upperBoundAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux)
{   
    setVariables(lowerBoundAux, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);
}


//QAGS method integral evaluator
double Integration1DimGSLQAGS::evaluate()
{   
    if ( fabs(lowerBound-upperBound)>minimumDistanceBetweenBounds )
    {
        //save original handler, turn off the error handler
        gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


        //allocate memory
        gsl_integration_workspace *workPtr = gsl_integration_workspace_alloc(workspaceLimitSize);
        
        int code = gsl_integration_qags(&F, lowerBound, upperBound, absolutePrecision, relativePrecision, workspaceLimitSize, workPtr, &result, &error);

        //handle possible integration error
        errorHandler(code, "QAGS");

        //clean memory
        gsl_integration_workspace_free(workPtr);

        
        //reset the error handler to previous state before switching off
        gsl_set_error_handler(old_error_handler);
    }
    else
    {
        result = 0.0;
    }

    return result;
}


//QAWC method constructor
Integration1DimGSLQAWC::Integration1DimGSLQAWC(double lowerBoundAux, double upperBoundAux, double singularityAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux)
{   
    setVariables(lowerBoundAux, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);

    //test if singularity is inside bounds
    if ( lowerBound<singularityAux && singularityAux<upperBound )
    {
        singularity = singularityAux;
    }
    else
    {
        cout << "Using QAWC method: singulariy is not inside the integration bounds!\n";
        abort();
    }
}


//QAWC method integral evaluator
double Integration1DimGSLQAWC::evaluate()
{   
    if ( fabs(lowerBound-upperBound)>minimumDistanceBetweenBounds )
    {
        //save original handler, turn off the error handler
        gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


        //allocate memory
        gsl_integration_workspace *workPtr = gsl_integration_workspace_alloc(workspaceLimitSize);
        
        int code = gsl_integration_qawc(&F, lowerBound, upperBound, singularity, absolutePrecision, relativePrecision, workspaceLimitSize, workPtr, &result, &error);

        //handle possible integration error
        errorHandler(code, "QAWC");

        //clean memory
        gsl_integration_workspace_free(workPtr);

        
        //reset the error handler to previous state before switching off
        gsl_set_error_handler(old_error_handler);
    }
    else
    {
        result = 0.0;
    }

    return result;
}


//QAGP method constructor
Integration1DimGSLQAGP::Integration1DimGSLQAGP(double lowerBoundAux, double upperBoundAux, vector<double> singularitiesAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux)
{   
    setVariables(lowerBoundAux, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);

    //check that all the provided singularities are inside the integration region
    for (int i = 0; i < int(singularitiesAux.size()); ++i)
    {
        if ( lowerBound<singularitiesAux[i] && singularitiesAux[i]<upperBound  )
        {
            singularities.push_back(singularitiesAux[i]);
        }
    }
    //QAGP reguires the singularities vector to contain the integration bounds
    singularities.push_back(lowerBound);
    singularities.push_back(upperBound);
    sort(singularities.begin(), singularities.end());
}


//QAGP method integral evaluator
double Integration1DimGSLQAGP::evaluate()
{   
    if ( fabs(lowerBound-upperBound)>minimumDistanceBetweenBounds )
    {
        //save original handler, turn off the error handler
        gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


        //allocate memory
        gsl_integration_workspace *workPtr = gsl_integration_workspace_alloc(workspaceLimitSize);
        
        int code = gsl_integration_qagp(&F, &singularities[0], int(singularities.size()), absolutePrecision, relativePrecision, workspaceLimitSize, workPtr,  &result, &error);

        //handle possible integration error
        errorHandler(code, "QAGP");

        //clean memory
        gsl_integration_workspace_free(workPtr);

        
        //reset the error handler to previous state before switching off
        gsl_set_error_handler(old_error_handler);
    }
    else
    {
        result = 0.0;
    }

    return result;
}


//CQUAD method constructor
Integration1DimGSLCQUAD::Integration1DimGSLCQUAD(double lowerBoundAux, double upperBoundAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux)
{   
    setVariables(lowerBoundAux, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);
}


//CQUAD method integral evaluator
double Integration1DimGSLCQUAD::evaluate()
{   
    if ( fabs(lowerBound-upperBound)>minimumDistanceBetweenBounds )
    {
        //save original handler, turn off the error handler
        gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


        //allocate memory
        gsl_integration_cquad_workspace *workPtr = gsl_integration_cquad_workspace_alloc(workspaceLimitSize);
        
        int code = gsl_integration_cquad(&F, lowerBound, upperBound, absolutePrecision, relativePrecision, workPtr, &result, &error, &numberFunctionEvaluations);

        //handle possible integration error
        errorHandler(code, "CQUAD");

        //clean memory
        gsl_integration_cquad_workspace_free(workPtr);

        
        //reset the error handler to previous state before switching off
        gsl_set_error_handler(old_error_handler);
    }
    else
    {
        result = 0.0;
    }

    return result;
}


//QAGI method constructor
Integration1DimGSLQAGI::Integration1DimGSLQAGI(GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux)
{   
    setVariables(0.0, 0.0, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);
}


//QAGI method integral evaluator
double Integration1DimGSLQAGI::evaluate()
{   
    //save original handler, turn off the error handler
    gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


    //allocate memory
    gsl_integration_workspace *workPtr = gsl_integration_workspace_alloc(workspaceLimitSize);
    
    int code = gsl_integration_qagi(&F, absolutePrecision, relativePrecision, workspaceLimitSize, workPtr, &result, &error);

    //handle possible integration error
    errorHandler(code, "QAGI");

    //clean memory
    gsl_integration_workspace_free(workPtr);

    
    //reset the error handler to previous state before switching off
    gsl_set_error_handler(old_error_handler);

    return result;
}


//QAGIU method constructor
Integration1DimGSLQAGIU::Integration1DimGSLQAGIU(double lowerBoundAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux)
{   
    setVariables(lowerBoundAux, 0.0, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);
}


//QAGIU method integral evaluator
double Integration1DimGSLQAGIU::evaluate()
{   
    //save original handler, turn off the error handler
    gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


    //allocate memory
    gsl_integration_workspace *workPtr = gsl_integration_workspace_alloc(workspaceLimitSize);
    
    int code = gsl_integration_qagiu(&F, lowerBound, absolutePrecision, relativePrecision, workspaceLimitSize, workPtr, &result, &error);

    //handle possible integration error
    errorHandler(code, "QAGIU");

    //clean memory
    gsl_integration_workspace_free(workPtr);

    
    //reset the error handler to previous state before switching off
    gsl_set_error_handler(old_error_handler);

    return result;
}


//QAGIL method constructor
Integration1DimGSLQAGIL::Integration1DimGSLQAGIL(double upperBoundAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux)
{   
    setVariables(0.0, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);
}


//QAGIL method integral evaluator
double Integration1DimGSLQAGIL::evaluate()
{   
    //save original handler, turn off the error handler
    gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


    //allocate memory
    gsl_integration_workspace *workPtr = gsl_integration_workspace_alloc(workspaceLimitSize);
    
    int code = gsl_integration_qagil(&F, upperBound, absolutePrecision, relativePrecision, workspaceLimitSize, workPtr, &result, &error);

    //handle possible integration error
    errorHandler(code, "QAGIL");

    //clean memory
    gsl_integration_workspace_free(workPtr);

    
    //reset the error handler to previous state before switching off
    gsl_set_error_handler(old_error_handler);

    return result;
}


//QAWS method constructor
Integration1DimGSLQAWS::Integration1DimGSLQAWS(double lowerBoundAux, double upperBoundAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux, double alphaAux, double betaAux, int muAux, int nuAux)
{   
    setVariables(lowerBoundAux, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);
    
    if ( alphaAux>-1.0 ){ alpha = alphaAux; }
    else{ cout << "Using QAWS method: problem with chosen alpha parameter! Remember, alpha>-1! \n "; abort(); }

    if ( betaAux>-1.0 ){ beta = betaAux; }
    else{ cout << "Using QAWS method: problem with chosen beta parameter! Remember, beta>-1! \n "; abort(); }

    if ( muAux==0 || muAux==1 ){ mu = muAux; }
    else{ cout << "Using QAWS method: problem with chosen mu parameter! Remember, mu=0 or mu=1! \n "; abort(); }

    if ( nuAux==0 || nuAux==1 ){ nu = nuAux; }
    else{ cout << "Using QAWS method: problem with chosen nu parameter! Remember, nu=0 or nu=1! \n "; abort(); }
}


//QAWS method integral evaluator
double Integration1DimGSLQAWS::evaluate()
{   
    if ( fabs(lowerBound-upperBound)>minimumDistanceBetweenBounds )
    {
        //save original handler, turn off the error handler
        gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off();


        //allocate memory
        gsl_integration_workspace *workPtr = gsl_integration_workspace_alloc(workspaceLimitSize);
        gsl_integration_qaws_table *table = gsl_integration_qaws_table_alloc(alpha, beta , mu, nu);
        
        int code = gsl_integration_qaws(&F, lowerBound, upperBound, table, absolutePrecision, relativePrecision, workspaceLimitSize, workPtr, &result, &error);

        //handle possible integration error
        errorHandler(code, "QAWS");

        //clean memory
        gsl_integration_qaws_table_free(table);
        gsl_integration_workspace_free(workPtr);

        
        //reset the error handler to previous state before switching off
        gsl_set_error_handler(old_error_handler);
    }
    else
    {
        result = 0.0;
    }
    
    return result;
}


double changeQAWCIntegrandToQAGSIntegrand(double x, void *parameters)
{   
    ChangeQAWCToQAGSParameters aux(parameters);
    double numerator = aux.evaluateNumerator(x);
    double singularity = aux.getSingularity();

    double function = numerator/( x - singularity );

    return function;
}


Integration1DimGSLQAWCQAGS::Integration1DimGSLQAWCQAGS(double lowerBoundAux, double upperBoundAux, double singularityAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux)
{   
    setVariables(lowerBoundAux, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);

    singularity = singularityAux;

    numeratorParameters = ChangeQAWCToQAGSParameters(integrandParameters, integrand, singularity);
}


Integration1DimGSLQAWCQAGS::Integration1DimGSLQAWCQAGS(double lowerBoundAux, double upperBoundAux, double singularityAux, GeneralIntegrandParameters* integrandParametersAux, double integrand(double, void*), double absolutePrecisionAux, double relativePrecisionAux, int workspaceLimitSizeAux , double minimumDistanceBetweenBoundAndSingularityAux)
{   
    setVariables(lowerBoundAux, upperBoundAux, integrandParametersAux, integrand, absolutePrecisionAux, relativePrecisionAux, workspaceLimitSizeAux);

    singularity = singularityAux;

    minimumDistanceBetweenBoundAndSingularity = minimumDistanceBetweenBoundAndSingularityAux;

    numeratorParameters = ChangeQAWCToQAGSParameters(integrandParameters, integrand, singularity);
}


bool Integration1DimGSLQAWCQAGS::isSingularityInsideTheIntegrationInterval()
{
    if ( fabs(lowerBound-upperBound)<minimumDistanceBetweenBounds )
    {
        //in this case we consider that lowerBound==upperBound!
        return false;
    }
    else
    {
        if ( lowerBound<upperBound )
        {
            if ( lowerBound<singularity && singularity<upperBound ){ return true; }
            else{ return false; }
        }
        else
        {
            if ( upperBound<singularity && singularity<lowerBound ){ return true; }
            else{ return false; }
        }
    }
}


double Integration1DimGSLQAWCQAGS::evaluate()
{   
    if ( isSingularityInsideTheIntegrationInterval()==true && 
         fabs(lowerBound-singularity)>minimumDistanceBetweenBoundAndSingularity && 
         fabs(singularity-upperBound)>minimumDistanceBetweenBoundAndSingularity )
    {   
        //evaluate the principal value of the integral
        Integration1DimGSLQAWC integralQAWC(lowerBound, upperBound, singularity, integrandParameters, F.function, absolutePrecision, relativePrecision, workspaceLimitSize);
        result = integralQAWC.evaluate();
    }
    else
    {   
        //evaluate the integral using the QAGS method
        Integration1DimGSLQAGS integralQAGS(lowerBound, upperBound, &numeratorParameters, changeQAWCIntegrandToQAGSIntegrand, absolutePrecision, relativePrecision, workspaceLimitSize);
        result = integralQAGS.evaluate();
    }

    return result;
}


double Integration1DimGSLQAWCQAGS::evaluateIntegration1DimNewtonCotes(int numberOfPartitions, NewtonCotesRule rule)
{   
    Integration1DimNewtonCotes trapezoidalSum(lowerBound, upperBound, numberOfPartitions, &numeratorParameters, changeQAWCIntegrandToQAGSIntegrand, rule);

    if ( isSingularityInsideTheIntegrationInterval()==true && 
         fabs(lowerBound-singularity)>minimumDistanceBetweenBoundAndSingularity && 
         fabs(singularity-upperBound)>minimumDistanceBetweenBoundAndSingularity )
    {   

        return trapezoidalSum.evaluateAvoidingSingularPoint(singularity);
    }
    else
    {   
        return trapezoidalSum.evaluate();
    }
}


////////////////////////////////////////////////////////////////////////////////////////


double integrandTest(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */
    double integrand = pow(x,2);

    return integrand;
}


double integrandTestCauchy(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */
    x = x; /* avoid unused parameter warning */

    double integrand = ( 1.0 );

    return integrand;
}


double integrandTestQAGP(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */

    double integrand = (1.0/sqrt(x))*log(x);

    return integrand;
}


double integrandTestQAGI(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */

    double integrand = exp( -pow(x,2) );

    return integrand;
}


double integrandRiemannCPV(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */

    double integrand = ( 1.0 )/( x-1.0 );

    return integrand;
}


double integrandTestQAWS(double x, void *parameters)
{   
    (void)(parameters); /* avoid unused parameter warning */

    double integrand = pow(x,2);

    return integrand;
}


void testIntegration1DimGSL()
{   
    cout << "Testing several GSL integration methods with different integrands.\n";
    cout << "All the integrals are normalized to 1.\n";

    double normalization = 0.0;

    TestIntegrandParameters aux1("integrandTest");
    Integration1DimGSLQNG integralQNG(-1.0, +2.0, &aux1, integrandTest, 1E-8, 1E-8);
    normalization = (1.0/3.0);
    double resultQNG = normalization*integralQNG.evaluate();
    cout << "resultQNG: " << resultQNG << "\n";


    Integration1DimGSLQAG integralQAG(-1.0, +2.0, &aux1, integrandTest, 1E-8, 1E-8, 1000, 1);
    normalization = (1.0/3.0);
    double resultQAG = normalization*integralQAG.evaluate();
    cout << "resultQAG: " << resultQAG << "\n";


    Integration1DimGSLQAGS integralQAGS(-1.0, +2.0, &aux1, integrandTest, 1E-8, 1E-8, 1000);
    normalization = (1.0/3.0);
    double resultQAGS = normalization*integralQAGS.evaluate();
    cout << "resultQAGS: " << resultQAGS << "\n";


    TestIntegrandParameters aux2("integrandTestCauchy");
    Integration1DimGSLQAWC integralQAWC(-1.0, +2.0, 1.0, &aux2, integrandTestCauchy, 1E-8, 1E-8, 1000);
    normalization = (-1.0/log(2.0));
    double resultQAWC = normalization*integralQAWC.evaluate();
    cout << "resultQAWC: " << resultQAWC << "\n";


    TestIntegrandParameters aux3("integrandTestQAGP");
    Integration1DimGSLQAGP integralQAGP(0, +1.0, {0.0}, &aux3, integrandTestQAGP, 1E-8, 1E-8, 1000);
    normalization = -0.25;
    double resultQAGP = normalization*integralQAGP.evaluate();
    cout << "resultQAGP: " << resultQAGP << "\n";

    Integration1DimGSLCQUAD integralCQUAD(0, +1.0, &aux3, integrandTestQAGP, 1E-8, 1E-8, 1000);
    normalization = -0.25;
    double resultCQUAD = normalization*integralCQUAD.evaluate();
    cout << "resultCQUAD: " << resultCQUAD << "\n";


    TestIntegrandParameters aux4("integrandTestQAGI");
    Integration1DimGSLQAGI integralQAGI(&aux4, integrandTestQAGI, 1E-8, 1E-8, 1000);
    normalization = ( sqrt(1.0/M_PI) );
    double resultQAGI = normalization*integralQAGI.evaluate();
    cout << "resultQAGI: " << resultQAGI << "\n";

    Integration1DimGSLQAGIU integralQAGIU(0, &aux4, integrandTestQAGI, 1E-8, 1E-8, 1000);
    normalization = ( sqrt(1.0/M_PI) );
    double resultQAGIU = normalization*2*integralQAGIU.evaluate();
    cout << "resultQAGIU: " << resultQAGIU << "\n";

    Integration1DimGSLQAGIL integralQAGIL(0, &aux4, integrandTestQAGI, 1E-8, 1E-8, 1000);
    normalization = ( sqrt(1.0/M_PI) );
    double resultQAGIL = normalization*2*integralQAGIL.evaluate();
    cout << "resultQAGIL: " << resultQAGIL << "\n";

    TestIntegrandParameters aux5("integrandTestQAWS");
    Integration1DimGSLQAWS integralQAWS(-1.0, +2.0, &aux1, integrandTestQAWS, 1E-8, 1E-8, 1000, -0.5, 0.0, 0, 0);
    normalization = ( 5.0/( 8.0*sqrt(3) ) );
    double resultQAWS = normalization*integralQAWS.evaluate();
    cout << "resultQAWS: " << resultQAWS << "\n";


    TestIntegrandParameters aux6("integrandTestQAWCQAGS");
    Integration1DimGSLQAWCQAGS integralQAWCQAGSIn(-5.0, 5.0, 4.0, &aux6, integrandTestCauchy, 1E-8, 1E-8, 1000);
    normalization = (-1.0/log(9.0));
    
    double resultQAWCQAGSIn = normalization*integralQAWCQAGSIn.evaluate();
    cout << "resultQAWCQAGSIn: " << resultQAWCQAGSIn << "\n";

    double compositeSumQAWCQAGSIn = normalization*integralQAWCQAGSIn.evaluateIntegration1DimNewtonCotes(10, alternativeCompositeSimpson);
    cout << "compositeSumQAWCQAGSIn: " << compositeSumQAWCQAGSIn << "\n";

    Integration1DimGSLQAWCQAGS integralQAWCQAGSOut(-5.0, 5.0, 8.0, &aux6, integrandTestCauchy, 1E-8, 1E-8, 1000);
    normalization = (-1.0/log(13.0/3.0));
    
    double resultQAWCQAGSOut = normalization*integralQAWCQAGSOut.evaluate();
    cout << "resultQAWCQAGSOut: " << resultQAWCQAGSOut << "\n";

    double compositeSumQAWCQAGSOut = normalization*integralQAWCQAGSOut.evaluateIntegration1DimNewtonCotes(10, alternativeCompositeSimpson);
    cout << "compositeSumQAWCQAGSOut: " << compositeSumQAWCQAGSOut << "\n";
}


void testIntegration1DimNewtonCotes()
{   
    cout << "Testing several Composite Trapezoidal Sum integration methods with different integrands.\n";
    cout << "All the integrals are normalized to 1.\n";

    double normalization = 0.0;

    TestIntegrandParameters aux1("integrandTest");
    Integration1DimNewtonCotes trapezoidalSum(-1.0, +2.0, 100, &aux1, integrandTest);
    normalization = (1.0/3.0);
    double resultTrapezoidalSum = normalization*trapezoidalSum.evaluate();
    cout << "resultTrapezoidalSum: " << resultTrapezoidalSum << "\n";

    TestIntegrandParameters aux2("integrandRiemannCPV");
    Integration1DimNewtonCotes trapezoidalSumCPV(-1.0, 2.0, 100, &aux2, integrandRiemannCPV, alternativeCompositeSimpson);
    normalization = (-1.0/log(2.0));
    double resultTrapezoidalSumCPV = normalization*trapezoidalSumCPV.evaluateAvoidingSingularPoint(1.0);
    cout << "resultTrapezoidalSumCPV: " << resultTrapezoidalSumCPV << "\n";
}



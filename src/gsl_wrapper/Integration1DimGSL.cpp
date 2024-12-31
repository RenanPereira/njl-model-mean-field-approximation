#include <iostream>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include "gsl_wrapper/Integration1DimGSL.h"

using namespace std;


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

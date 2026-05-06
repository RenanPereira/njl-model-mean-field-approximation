#include <iostream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff_klev_recipe.h"
#include "njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff.h"
#include "physics_utils/distribution_functions.h"

using namespace std;


//This function is called if an error is found during a PRO integration.
//If the integral_ID is present in this function, the parameters that are being used in the function
//that is being integrated are presented in the console for debugging purposes.
double integral_inspector(string integral_ID, void* parameters, int code)
{   
    //Print the type of error
    if      (code==GSL_EMAXITER){ printf ("%s\n", gsl_strerror (GSL_EMAXITER)); }
    else if (code==GSL_EDIVERGE){ printf ("%s\n", gsl_strerror (GSL_EDIVERGE)); }
    else if (code==GSL_ESING)   { printf ("%s\n", gsl_strerror (GSL_ESING)); }
    else if (code==GSL_EROUND)  { printf ("%s\n", gsl_strerror (GSL_EROUND)); }

    //Depending on the integral ID, return the parameters that enter the function
    if( integral_ID=="IntB0RehbergKlevanskypkfinTfinRe" || integral_ID=="IntB0RehbergKlevanskymkfinTfinRe" )
    {
        double Lambda = ((struct f1_loop_parameters *)(parameters))->cutoff;
        double T = ((struct f1_loop_parameters *)(parameters))->temperature;
        double Cp = ((struct f1_loop_parameters *)(parameters))->eff_chem_pot_1;
        double Mi=((struct f1_loop_parameters *)(parameters))->eff_mass_quark_1;
        double Mj=((struct f1_loop_parameters *)(parameters))->eff_mass_quark_2;
        double k=((struct f1_loop_parameters *)(parameters))->external_momentum;
        double lambdax=((struct f1_loop_parameters *)(parameters))->omega;
        printf("Parameters:(lambdax, Mi, Mj, k, Lambda, T, Cp) \n");
        printf("%.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f \n", lambdax, Mi, Mj, k, Lambda,T,Cp);
    }

    //In the end, abort the code
    abort();

    return 0;
}

//QAGS ADAPTIVE INTEGRATION WITH SINGULARITIES
//This function applies the Gauss-Kronrod 21-point integration rule adaptively. The results are extrapolated using the
//epsilon-algorithm, which accelerates the convergence of the integral in the presence of discontinuities and integrable singularities.
//The subintervals and their results are stored in the memory provided by workspace. The maximum number of subintervals is 
//given by limit, which may not exceed the allocated size of the workspace.
double integrate_QAGS_PRO(string integral_ID, double x_init, double x_final, void* params, double placeholder_f (double, void*), int workspace_size, double abs_precision, double rel_precision)
{   
    gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off(); //save original handler, turn off the error handler

    const size_t integration_workspace_size = workspace_size;
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(integration_workspace_size);

    gsl_function F;
    F.function = placeholder_f;
    F.params = params;

    double result, error;

    const double xlow = x_init;
    const double xhigh = x_final;
    const double epsabs = abs_precision;
    const double epsrel = rel_precision;
    
    int code = gsl_integration_qags(&F, xlow, xhigh, epsabs, epsrel, integration_workspace_size, work_ptr, &result, &error);
    
    if( code!=0 )
    {   
        cout << "Problem in the integration using QAGS: " << integral_ID << "\n";
        integral_inspector(integral_ID, params, code);
    }

    gsl_integration_workspace_free(work_ptr);

    gsl_set_error_handler(old_error_handler); //reset the error handler 
    
    return result;
}

//QAGP
//This function applies the adaptive integration algorithm QAGS taking account of the user-supplied locations of singular points. 
//The array pts of length npts should contain the endpoints of the integration ranges deﬁned by the integration region and locations
//of the singularities.
double integrate_QAGP_PRO(string integral_ID, double x_init, vector<double> singular_points, double x_final, void* params, double placeholder_f (double, void*), int workspace_size, double abs_precision, double rel_precision)
{
    gsl_error_handler_t *old_error_handler = gsl_set_error_handler_off(); //save original handler, turn off the error handler

    const size_t integration_workspace_size = workspace_size;
    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(integration_workspace_size);

    gsl_function F;
    F.function = placeholder_f;
    F.params = params;

    double result, error;

    const double xlow = x_init;
    const double xhigh = x_final;
    const vector<double> xsing = singular_points;
    const double epsabs = abs_precision;
    const double epsrel = rel_precision;

    vector<double> x;

    x.push_back(xlow);
    for (int i=0; i< int(singular_points.size()); i++)
    {
        x.push_back(singular_points[i]);
    }
    x.push_back(xhigh);
    sort(x.begin(), x.end());

    const size_t xsize = 2 + singular_points.size();

    int code = gsl_integration_qagp (&F, &x[0], xsize, epsabs, epsrel, 1000, work_ptr,  &result, &error);

    if( code!=0 )
    {   
        cout << "Problem in the integration using QAGP: " << integral_ID << "\n";
        integral_inspector(integral_ID, params, code);
    }

    gsl_integration_workspace_free(work_ptr);

    gsl_set_error_handler(old_error_handler); //reset the error handler 
    
    return result;
}


////////////////////////////////////////////////////////////


//We introduce the chemical potential a la Klevansky: wn -> wn-i*mu
//For the imaginary part, we use the Klevansky shift. The pole is shifted as follows: eta*omega -> eta*omega - i*epsilon


//Parameters necessary in this file
double f1_Klev_precision = 1E-6;
int int_Klev_workspace = 2000;
double ZERO_Klev=1E-15;


//////////////////////////////////////////////////////////////////
//Integral B Klevansky
//////////////////////////////////////////////////////////////////
double Intfermidist(double T,double EnerA,double EnerB){
    double aux=EnerB-EnerA-T*log( 1. + exp(EnerB/T))+T*log( 1. + exp(EnerA/T));
    return aux;
}

double IntfermidistT0(double EnerA,double EnerB){
    double aux= EnerB-EnerA+
                heavisideTheta(EnerA)*heavisideTheta(EnerA-EnerB)*(EnerA - EnerB*heavisideTheta(EnerB))+
                (-EnerB +EnerA*heavisideTheta(EnerA))*heavisideTheta(EnerB)*heavisideTheta(-EnerA+EnerB);
    return aux;
}
//////////////////////
double doublesign(double x){
    if (x > 0) return +1.;
    if (x < 0) return -1.;
    return 0.;
}
////////////////////////////////////////////
double funEner1lambdax0(double k, double Mi, double Mj){
    double aux;
    aux=sqrt(pow(pow(k,2) + pow(Mi,2),2) + 2*(k - Mi)*(k + Mi)*pow(Mj,2) + pow(Mj,4))/(2.*k);
    return aux;
    }

double funEner1lambdaxnot0(double lambdax, double k, double Mi, double Mj){
    double aux;
    aux=(lambdax*(-pow(k,2) + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2)))/(2.*(pow(k,2) - pow(lambdax,2))) + 
   (k*sqrt(4*(pow(k,2) - pow(lambdax,2))*pow(Mi,2) + pow(-pow(k,2) + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2),2)))/(2.*(pow(k,2) - pow(lambdax,2)));
    return aux;
    }

double funEner1( double lambdax, double k, double Mi, double Mj){
    double aux;
    if(fabs(lambdax)<=ZERO_Klev){
        aux=funEner1lambdax0(k,Mi,Mj);
    }else{
        aux=funEner1lambdaxnot0(lambdax,k,Mi,Mj);    
    }
    return aux;
    }
//////////////////////
double funEner2lambdax0(double k, double Mi, double Mj){
    double aux;
    aux=-sqrt(pow(pow(k,2) + pow(Mi,2),2) + 2*(k - Mi)*(k + Mi)*pow(Mj,2) + pow(Mj,4))/(2.*k);
    return aux;
    }

double funEner2lambdaxnot0(double lambdax, double k, double Mi, double Mj){
    double aux;
    aux=(lambdax*(-pow(k,2) + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2)))/(2.*(pow(k,2) - pow(lambdax,2))) - 
   (k*sqrt(4*(pow(k,2) - pow(lambdax,2))*pow(Mi,2) + pow(-pow(k,2) + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2),2)))/(2.*(pow(k,2) - pow(lambdax,2)));
    return aux;
    }

double funEner2( double lambdax, double k, double Mi, double Mj){
    double aux;
    if(fabs(lambdax)<=ZERO_Klev){
        aux=funEner2lambdax0(k,Mi,Mj);
    }else{
        aux=funEner2lambdaxnot0(lambdax,k,Mi,Mj);    
    }
    return aux;
}
//////////////////////
double funEnerlambdaxkequal( double k, double Mi, double Mj){
    double aux;
    aux=-k*pow(Mi,2)/(pow(Mi,2)-pow(Mj,2))-(pow(Mi,2)-pow(Mj,2))/(4*k);
    return aux;
}
////////////////////////////////////////////
// double IntegrandoReIntB0RehbergKlevanskypEner(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,+Ener-Cp)
//           *log(pow(+1. + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2)/
//                pow(-1. + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2));
// return aux;
// }
// double IntegrandoReIntB0RehbergKlevanskymEner(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,-Ener-Cp)
//           *log(pow(1 + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2)/
//                pow(-1 + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2));
// return aux;
// }

//////////
// double IntegrandoReIntB0RehbergKlevanskypEner(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,+Ener-Cp)
//           *(+log(pow(+1. + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2))
//             -log(pow(-1. + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2))
//             );
// return aux;
// }

// double IntegrandoReIntB0RehbergKlevanskymEner(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,-Ener-Cp)
//           *(
//             +log(pow(+1. + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2) )/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2))
//             -log(pow(-1. + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2))
//             );
//     return aux;
// }

//////////
double IntegrandoReIntB0RehbergKlevanskypEner(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
{
double aux=fermiDistribution(T,+Ener-Cp)
          *(    
                +log(pow(+2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2)),2))
                -log(pow(-2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2)),2))
            );
return aux;
}

double IntegrandoReIntB0RehbergKlevanskymEner(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
{
double aux=fermiDistribution(T,-Ener-Cp)
          *(    
                +log(pow(+2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2)),2))
                -log(pow(-2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2)),2))
            );

return aux;
}

// double IntegrandoReIntB0RehbergKlevanskypEnerA(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,+Ener-Cp)
//           *(+log(pow(+2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2)),2)));

//     return aux;
// }

// double IntegrandoReIntB0RehbergKlevanskypEnerB(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,+Ener-Cp)
//           *(-log(pow(-2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2)),2)));
//     return aux;
// }

// double IntegrandoReIntB0RehbergKlevanskymEnerA(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
//     double aux=fermiDistribution(T,-Ener-Cp)
//           *(+log(pow(+2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2)),2)));
//     return aux;
// }

// double IntegrandoReIntB0RehbergKlevanskymEnerB(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,-Ener-Cp)
//           *(-log(pow(-2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2)),2)));
//     return aux;
// }

// double IntegrandoReIntB0RehbergKlevanskypEnerA(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,+Ener-Cp)
//           *(+log(pow(+1. + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2)));

//     return aux;
// }

// double IntegrandoReIntB0RehbergKlevanskypEnerB(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,+Ener-Cp)
//           *(-log(pow(-1. + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2)));
//     return aux;
// }

// double IntegrandoReIntB0RehbergKlevanskymEnerA(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
//     double aux=fermiDistribution(T,-Ener-Cp)
//           *(+log(pow(+1. + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2)));
//     return aux;
// }

// double IntegrandoReIntB0RehbergKlevanskymEnerB(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
// {
// double aux=fermiDistribution(T,-Ener-Cp)
//           *(-log(pow(-1. + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2)));
//     return aux;
// }
double IntegrandoReIntB0RehbergKlevanskypEnerA(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
{
double aux=fermiDistribution(T,+Ener-Cp)
          *(+log(pow(+2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2)),2)));

    return aux;
}

double IntegrandoReIntB0RehbergKlevanskypEnerB(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
{
double aux=fermiDistribution(T,+Ener-Cp)
          *(-log(pow(-2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2)),2)));
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskymEnerA(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
{
    double aux=fermiDistribution(T,-Ener-Cp)
          *(+log(pow(+2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2)),2)));
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskymEnerB(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double Ener)
{
double aux=fermiDistribution(T,-Ener-Cp)
          *(-log(pow(-2.*k*sqrt(pow(Ener,2) - pow(Mi,2)) + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2)),2)));
    return aux;
}
double IntegrandoReIntB0RehbergKlevanskypp(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double p)
{
//      double aux=fermiDistribution(T,+sqrt(pow(Mi,2)+pow(p,2))-Cp)
//           *p/sqrt(pow(Mi,2)+pow(p,2))*(
//             +log(pow(+2.*k*p + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*sqrt(pow(Mi,2)+pow(p,2))*lambdax + pow(lambdax,2)),2))
//             -log(pow(-2.*k*p + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*sqrt(pow(Mi,2)+pow(p,2))*lambdax + pow(lambdax,2)),2))
//             );
    double aux=p/sqrt(pow(Mi,2)+pow(p,2))*IntegrandoReIntB0RehbergKlevanskypEner(lambdax,Mi,Mj,k,Lambda,T,Cp,sqrt(pow(Mi,2)+pow(p,2)));
return aux;
}

double IntegrandoReIntB0RehbergKlevanskymp(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double p)
{
    // double aux=fermiDistribution(T,-sqrt(pow(Mi,2)+pow(p,2))-Cp)
    //       *p/sqrt(pow(Mi,2)+pow(p,2))*(
    //         +log(pow(+2.*k*p + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*sqrt(pow(Mi,2)+pow(p,2))*lambdax + pow(lambdax,2)),2))
    //         -log(pow(-2.*k*p + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*sqrt(pow(Mi,2)+pow(p,2))*lambdax + pow(lambdax,2)),2))
    //         );
    double aux=p/sqrt(pow(Mi,2)+pow(p,2))*IntegrandoReIntB0RehbergKlevanskymEner(lambdax,Mi,Mj,k,Lambda,T,Cp,sqrt(pow(Mi,2)+pow(p,2)));

return aux;
}

double IntegrandoReIntB0RehbergKlevanskyppA(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double p)
{
    double aux=p/sqrt(pow(Mi,2)+pow(p,2))*IntegrandoReIntB0RehbergKlevanskypEnerA(lambdax,Mi,Mj,k,Lambda,T,Cp,sqrt(pow(Mi,2)+pow(p,2)));
return aux;
}

double IntegrandoReIntB0RehbergKlevanskyppB(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double p)
{
    double aux=p/sqrt(pow(Mi,2)+pow(p,2))*IntegrandoReIntB0RehbergKlevanskypEnerB(lambdax,Mi,Mj,k,Lambda,T,Cp,sqrt(pow(Mi,2)+pow(p,2)));
return aux;
}

double IntegrandoReIntB0RehbergKlevanskympA(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double p)
{
    double aux=p/sqrt(pow(Mi,2)+pow(p,2))*IntegrandoReIntB0RehbergKlevanskymEnerA(lambdax,Mi,Mj,k,Lambda,T,Cp,sqrt(pow(Mi,2)+pow(p,2)));

return aux;
}

double IntegrandoReIntB0RehbergKlevanskympB(double lambdax, double Mi, double Mj,double k,double Lambda,double T, double Cp,double p)
{
    double aux=p/sqrt(pow(Mi,2)+pow(p,2))*IntegrandoReIntB0RehbergKlevanskymEnerB(lambdax,Mi,Mj,k,Lambda,T,Cp,sqrt(pow(Mi,2)+pow(p,2)));

return aux;
}

double IntegrandoReIntB0RehbergKlevanskypEner1D(double Ener, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskypEner(lambdax, Mi, Mj, k, Lambda, T, Cp, Ener);
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskymEner1D(double Ener, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskymEner(lambdax, Mi, Mj, k, Lambda, T, Cp, Ener);
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskypp1D(double p, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskypp(lambdax, Mi, Mj, k, Lambda, T, Cp, p);
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskymp1D(double p, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskymp(lambdax, Mi, Mj, k, Lambda, T, Cp, p);
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskypEnerA1D(double Ener, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskypEnerA(lambdax, Mi, Mj, k, Lambda, T, Cp, Ener);
    return aux;
}
double IntegrandoReIntB0RehbergKlevanskypEnerB1D(double Ener, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskypEnerB(lambdax, Mi, Mj, k, Lambda, T, Cp, Ener);
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskymEnerA1D(double Ener, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskymEnerA(lambdax, Mi, Mj, k, Lambda, T, Cp, Ener);
    return aux;
}
double IntegrandoReIntB0RehbergKlevanskymEnerB1D(double Ener, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskymEnerB(lambdax, Mi, Mj, k, Lambda, T, Cp, Ener);
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskyppA1D(double p, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskyppA(lambdax, Mi, Mj, k, Lambda, T, Cp, p);
    return aux;
}
double IntegrandoReIntB0RehbergKlevanskyppB1D(double p, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskyppB(lambdax, Mi, Mj, k, Lambda, T, Cp, p);
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskympA1D(double p, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskympA(lambdax, Mi, Mj, k, Lambda, T, Cp, p);
    return aux;
}
double IntegrandoReIntB0RehbergKlevanskympB1D(double p, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskympB(lambdax, Mi, Mj, k, Lambda, T, Cp, p);
    return aux;
}


double IntegrandoReIntB0RehbergKlevanskypT0Ener(double lambdax, double Mi, double Mj,double k,double Lambda, double Cp,double Ener)
{
double aux=(1.-heavisideTheta(+Ener-Cp))
          *log(pow(1 + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2)/
               pow(-1 + (-pow(k,2) + pow(Mi,2) - pow(Mj,2) + 2*Ener*lambdax + pow(lambdax,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2))
;
return aux;
}

double IntegrandoReIntB0RehbergKlevanskymT0Ener(double lambdax, double Mi, double Mj,double k,double Lambda, double Cp,double Ener)
{
double aux=(1.-heavisideTheta(-Ener-Cp))
          *log(pow(1 + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2)/
               pow(-1 + (-pow(k,2) + 2*Ener*lambdax + pow(lambdax,2) + pow(Mi,2) - pow(Mj,2))/(2.*k*sqrt(pow(Ener,2) - pow(Mi,2))),2));

return aux;
}

double IntegrandoReIntB0RehbergKlevanskypT0Enerlambdax0(double Mi, double Mj,double k,double Lambda, double Cp,double Ener)
{
double aux=(1.-heavisideTheta(+Ener-Cp))
          *log(pow(pow(k,2) - pow(Mi,2) - 2*k*sqrt((Ener - Mi)*(Ener + Mi)) + pow(Mj,2),2)/
               pow(pow(k,2) - pow(Mi,2) + 2*k*sqrt((Ener - Mi)*(Ener + Mi)) + pow(Mj,2),2));
return aux;
}

double IntegrandoReIntB0RehbergKlevanskymT0Enerlambdax0(double Mi, double Mj,double k,double Lambda, double Cp,double Ener)
{
double aux=(1.-heavisideTheta(-Ener-Cp))
          *log(pow(pow(k,2) - pow(Mi,2) - 2*k*sqrt((Ener - Mi)*(Ener + Mi)) + pow(Mj,2),2)/
               pow(pow(k,2) - pow(Mi,2) + 2*k*sqrt((Ener - Mi)*(Ener + Mi)) + pow(Mj,2),2))
          ;
return aux;
}

double IntegrandoReIntB0RehbergKlevanskypT0Ener2(double lambdax, double Mi, double Mj,double k,double Lambda, double Cp,double Ener)
{
double aux;
if(fabs(lambdax)<=ZERO_Klev){
    aux=IntegrandoReIntB0RehbergKlevanskypT0Enerlambdax0(Mi,Mj,k,Lambda,Cp,Ener);
    }
else{
     aux=IntegrandoReIntB0RehbergKlevanskypT0Ener(lambdax,Mi,Mj,k,Lambda,Cp,Ener);
}
return aux;
}

double IntegrandoReIntB0RehbergKlevanskymT0Ener2(double lambdax, double Mi, double Mj,double k,double Lambda, double Cp,double Ener)
{
double aux;
if(fabs(lambdax)<=ZERO_Klev){
    aux=IntegrandoReIntB0RehbergKlevanskymT0Enerlambdax0(Mi,Mj,k,Lambda,Cp,Ener);
    }
else{
     aux=IntegrandoReIntB0RehbergKlevanskymT0Ener(lambdax,Mi,Mj,k,Lambda,Cp,Ener);
}
return aux;
}

double IntegrandoReIntB0RehbergKlevanskypT0Ener21D(double Ener, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskypT0Ener2(lambdax, Mi, Mj, k, Lambda, Cp, Ener);
    return aux;
}

double IntegrandoReIntB0RehbergKlevanskymT0Ener21D(double Ener, void* params)
{
    double Mi=((struct f1_loop_parameters*)(params))->eff_mass_quark_1;
    double Mj=((struct f1_loop_parameters*)(params))->eff_mass_quark_2;
    double Lambda=((struct f1_loop_parameters*)(params))->cutoff;
    double T=((struct f1_loop_parameters*)(params))->temperature;
    double Cp=((struct f1_loop_parameters*)(params))->eff_chem_pot_1;
    double k=((struct f1_loop_parameters*)(params))->external_momentum;
    double lambdax=((struct f1_loop_parameters*)(params))->omega;
    double aux=IntegrandoReIntB0RehbergKlevanskymT0Ener2(lambdax, Mi, Mj, k, Lambda, Cp, Ener);
    return aux;
}

//////////////////////
vector<double> polos(double lambdax, double Mi, double Mj, double k, double Lambda){
    vector<double> auxpolos={};

    if(fabs(lambdax-(+k))>ZERO_Klev && fabs(lambdax-(-k))>ZERO_Klev)
    {
        if(funEner1(lambdax, k, Mi, Mj)>=Mi && funEner1(lambdax, k, Mi, Mj)<=sqrt(pow(Mi,2)+pow(Lambda,2))){
            auxpolos.push_back(funEner1(lambdax, k, Mi, Mj));}
        if(funEner2(lambdax, k, Mi, Mj)>=Mi && funEner2(lambdax, k, Mi, Mj)<=sqrt(pow(Mi,2)+pow(Lambda,2))){
            auxpolos.push_back(funEner2(lambdax, k, Mi, Mj));}}
    else
    {
        if (fabs(Mi-Mj)>ZERO_Klev)
        {
            if(fabs(lambdax-k)<=ZERO_Klev)
            {
                if(+funEnerlambdaxkequal( k, Mi, Mj)>=Mi && +funEnerlambdaxkequal( k, Mi, Mj)<=sqrt(pow(Mi,2)+pow(Lambda,2)))
                {
                    auxpolos.push_back(+funEnerlambdaxkequal( k, Mi, Mj));
                }
            }else
            {
                if(-funEnerlambdaxkequal( k, Mi, Mj)>=Mi && -funEnerlambdaxkequal( k, Mi, Mj)<=sqrt(pow(Mi,2)+pow(Lambda,2)))
                {
                    auxpolos.push_back(-funEnerlambdaxkequal( k, Mi, Mj));
                }
            }
        }
    }    
    return auxpolos;
}
///////////For k=0 results we will use the definitions coming from f1_loop_function so these are just dummy functions
//////////////////////
double IntB0RehbergKlevanskypk0T0Re(double lambdax,double Mi, double Mj,double Lambda, double Cp){
    double aux=0;
    return aux;
}

double IntB0RehbergKlevanskypk0T0Im(double lambdax,double Mi, double Mj,double Lambda, double Cp){
    double aux=0;
    return aux;
}

double IntB0RehbergKlevanskymk0T0Re(double lambdax,double Mi, double Mj,double Lambda, double Cp){
    double aux=0;
    return aux;
}

double IntB0RehbergKlevanskymk0T0Im(double lambdax,double Mi, double Mj,double Lambda, double Cp){
    double aux=0;
    return aux;
}

double IntB0RehbergKlevanskypk0TfinRe(double lambdax,double Mi, double Mj,double Lambda, double T,double Cp){
    double aux=0;
    return aux;
}

double IntB0RehbergKlevanskypk0TfinIm(double lambdax,double Mi, double Mj,double Lambda, double T,double Cp){
    double aux=0;
    return aux;
}

double IntB0RehbergKlevanskymk0TfinRe(double lambdax,double Mi, double Mj,double Lambda, double T,double Cp){
    double aux=0;
    return aux;
}

double IntB0RehbergKlevanskymk0TfinIm(double lambdax,double Mi, double Mj,double Lambda, double T,double Cp){
    double aux=0;
    return aux;
}
//////////////////////
double IntB0RehbergKlevanskypkfinT0Re(double lambdax,double Mi,double Mj,double k,double Lambda,double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = 0;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;
    vector<double> auxpolos={};
    
    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);

    if (auxpolos.size()==0){
        aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskypkfinT0Re", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, & IntegrandoReIntB0RehbergKlevanskypT0Ener21D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }else{
        aux=1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskypkfinT0Re", Mi, auxpolos, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskypT0Ener21D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }
    return aux;
}

double IntB0RehbergKlevanskymkfinT0Re(double lambdax,double Mi,double Mj,double k,double Lambda,double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = 0;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;
    vector<double> auxpolos={};

    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);

    if (auxpolos.size()==0)
    {
        aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskymkfinT0Re", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, & IntegrandoReIntB0RehbergKlevanskymT0Ener21D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }else
    {
        aux=1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskymkfinT0Re", Mi, auxpolos, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskymT0Ener21D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }
    return aux;
}

double IntB0RehbergKlevanskypkfinT0Im(double lambdax,double Mi,double Mj,double k,double Lambda,double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = 0;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;
    vector<double> auxpolos={};
    
    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);
    
    if (auxpolos.size()==0){
        aux=0;
        }else{
            if(auxpolos.size()==1){
                if(fabs(lambdax-(+k))>ZERO_Klev && fabs(lambdax-(-k))>ZERO_Klev)
                {
                    if(funEner1(lambdax, k, Mi, Mj)>=Mi && funEner1(lambdax, k, Mi, Mj)<=sqrt(pow(Mi,2)+pow(Lambda,2))){
                        aux=M_PI/k*doublesign(lambdax)*IntfermidistT0(funEner1(lambdax, k, Mi, Mj)-Cp,sqrt(pow(Mi,2)+pow(Lambda,2))-Cp);
                    }else{
                        aux=M_PI/k*doublesign(lambdax)*IntfermidistT0(funEner2(lambdax, k, Mi, Mj)-Cp,sqrt(pow(Mi,2)+pow(Lambda,2))-Cp);
                    }
                }else{
                    if(fabs(Mi-Mj)>ZERO_Klev)
                    {
                        if(fabs(lambdax-(+k))<=ZERO_Klev)
                        {
                            aux=M_PI/k*doublesign(lambdax)*IntfermidistT0(+funEnerlambdaxkequal( k, Mi, Mj)-Cp,sqrt(pow(Mi,2)+pow(Lambda,2))-Cp);
                        }else{
                            aux=M_PI/k*doublesign(lambdax)*IntfermidistT0(-funEnerlambdaxkequal( k, Mi, Mj)-Cp,sqrt(pow(Mi,2)+pow(Lambda,2))-Cp);
                        }
                    }
                }
            }else{
                if(auxpolos.size()==2)
            {
                aux=M_PI/k*doublesign(lambdax)*
                    IntfermidistT0(
                    min(funEner1(lambdax, k, Mi, Mj),funEner2(lambdax, k, Mi, Mj))-Cp,
                    max(funEner1(lambdax, k, Mi, Mj),funEner2(lambdax, k, Mi, Mj))-Cp);
            }
        }
     }                   
    return aux;
}

double IntB0RehbergKlevanskymkfinT0Im(double lambdax,double Mi,double Mj,double k,double Lambda,double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = 0;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;
    vector<double> auxpolos={};

    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);
    if (auxpolos.size()==0){
        aux=0;
        }else{
            if(auxpolos.size()==1){
                if(fabs(lambdax-(+k))>ZERO_Klev && fabs(lambdax-(-k))>ZERO_Klev )
                {
                    if(funEner1(lambdax, k, Mi, Mj)>=Mi && funEner1(lambdax, k, Mi, Mj)<=sqrt(pow(Mi,2)+pow(Lambda,2))){
                        aux=-M_PI/k*doublesign(lambdax)*IntfermidistT0(-(funEner1(lambdax, k, Mi, Mj)+Cp),-(sqrt(pow(Mi,2)+pow(Lambda,2))+Cp));
                    }else{
                        aux=-M_PI/k*doublesign(lambdax)*IntfermidistT0(-(funEner2(lambdax, k, Mi, Mj)+Cp),-(sqrt(pow(Mi,2)+pow(Lambda,2))+Cp));                    }
                }else{
                    if(fabs(Mi-Mj)>ZERO_Klev)
                    {
                        if(fabs(lambdax-(+k))<=ZERO_Klev)
                        {
                            aux=-M_PI/k*doublesign(lambdax)*IntfermidistT0(-(funEnerlambdaxkequal(k, Mi, Mj)+Cp),-(sqrt(pow(Mi,2)+pow(Lambda,2))+Cp));
                        }else{
                            aux=-M_PI/k*doublesign(lambdax)*IntfermidistT0(-(-funEnerlambdaxkequal(k, Mi, Mj)+Cp),-(sqrt(pow(Mi,2)+pow(Lambda,2))+Cp));
                        }
                    }
                }
            }else{
                if(auxpolos.size()==2)
            {
                aux=M_PI/k*doublesign(lambdax)*
                    IntfermidistT0(
                    -(min(funEner1(lambdax, k, Mi, Mj),funEner2(lambdax, k, Mi, Mj))+Cp),
                    -(max(funEner1(lambdax, k, Mi, Mj),funEner2(lambdax, k, Mi, Mj))+Cp));
            }
        }
     }                   
    return aux;
}

////////////////////
double IntB0RehbergKlevanskypkfinTfinReAlt(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = T;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;
    vector<double> auxpolos={};

    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);

    if (auxpolos.size()==0){
        aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskypkfinTfinRe", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskypEner1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }else{
        aux=1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskypkfinTfinRe", Mi, auxpolos, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskypEner1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }   
    return aux;
}

double IntB0RehbergKlevanskymkfinTfinReAlt(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = T;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;
    vector<double> auxpolos={};

    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);

    if (auxpolos.size()==0){
        aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskymkfinTfinRe", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, & IntegrandoReIntB0RehbergKlevanskymEner1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }else{
        aux=1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskymkfinTfinRe", Mi, auxpolos, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskymEner1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }   
    return aux;
}


// double IntB0RehbergKlevanskypkfinTfinReAlt(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
//     double aux=0;
//     struct f1_loop_parameters auxparams;
//     auxparams.temperature = T;
//     auxparams.eff_chem_pot_1 = Cp;
//     auxparams.cutoff = Lambda;
//     auxparams.eff_mass_quark_1 = Mi;
//     auxparams.eff_mass_quark_2 = Mj;
//     auxparams.omega = lambdax;
//     auxparams.external_momentum=k;
    
//     vector<double> auxpolos={};
//     auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);

//     vector<double> auxpolosp={};
//     for (int i=0;i<auxpolos.size();i++){
//         auxpolosp.push_back( sqrt(pow(auxpolos[i],2)-pow(Mi,2)));
//     }

//     if (auxpolos.size()==0){
//         aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskypkfinTfinRe", 0, Lambda, &auxparams, &IntegrandoReIntB0RehbergKlevanskypp1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     }else{
//         aux=1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskypkfinTfinRe", 0, auxpolosp, Lambda, &auxparams, &IntegrandoReIntB0RehbergKlevanskypp1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     }   
//     return aux;
// }

// double IntB0RehbergKlevanskymkfinTfinReAlt(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
//     double aux=0;
//     struct f1_loop_parameters auxparams;
//     auxparams.temperature = T;
//     auxparams.eff_chem_pot_1 = Cp;
//     auxparams.cutoff = Lambda;
//     auxparams.eff_mass_quark_1 = Mi;
//     auxparams.eff_mass_quark_2 = Mj;
//     auxparams.omega = lambdax;
//     auxparams.external_momentum=k;
//     vector<double> auxpolos={};

//     auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);
    
//     vector<double> auxpolosp={};
//     for (int i=0;i<auxpolos.size();i++){
//         auxpolosp.push_back( sqrt(pow(auxpolos[i],2)-pow(Mi,2)));
//     }

//     if (auxpolos.size()==0){
//         aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskymkfinTfinRe", 0, Lambda, &auxparams, & IntegrandoReIntB0RehbergKlevanskymp1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     }else{
//         aux=1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskymkfinTfinRe", 0, auxpolosp, Lambda, &auxparams, &IntegrandoReIntB0RehbergKlevanskymp1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     }   
//     return aux;
// }

// double IntB0RehbergKlevanskypkfinTfinRe(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
//     double aux=0;
//     struct f1_loop_parameters auxparams;
//     auxparams.temperature = T;
//     auxparams.eff_chem_pot_1 = Cp;
//     auxparams.cutoff = Lambda;
//     auxparams.eff_mass_quark_1 = Mi;
//     auxparams.eff_mass_quark_2 = Mj;
//     auxparams.omega = lambdax;
//     auxparams.external_momentum=k;
//     vector<double> auxpolos={};

//     auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);
    
//     aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskypkfinTfinRe", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskypEnerB1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     cout<<"teste "<<aux<<"\n";
//     cout<<"num de polos: "<< auxpolos.size()<<"\n";
//     if (auxpolos.size()>0){
//         cout<<"polo: "<< auxpolos[0]<< "  " <<IntegrandoReIntB0RehbergKlevanskypEnerA1D(auxpolos[0], &auxparams)<<"\n";
//         }

//     if (auxpolos.size()==0){
//         // aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskypkfinTfinRe", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskypEnerB1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//         aux=aux+1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskypkfinTfinRe", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskypEnerA1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     }else{
//         // aux=1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskypkfinTfinRe", Mi, auxpolos, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskypEnerB1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//         aux=aux+1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskypkfinTfinRe", Mi, auxpolos, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskypEnerA1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     }   
//     return aux;
// }

// double IntB0RehbergKlevanskymkfinTfinRe(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
//     double aux=0;
//     struct f1_loop_parameters auxparams;
//     auxparams.temperature = T;
//     auxparams.eff_chem_pot_1 = Cp;
//     auxparams.cutoff = Lambda;
//     auxparams.eff_mass_quark_1 = Mi;
//     auxparams.eff_mass_quark_2 = Mj;
//     auxparams.omega = lambdax;
//     auxparams.external_momentum=k;
//     vector<double> auxpolos={};

//     auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);
    
//     aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskymkfinTfinRe", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskymEnerB1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     if (auxpolos.size()==0){
//         // aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskymkfinTfinRe", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskymEnerB1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//         aux=aux+1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskymkfinTfinRe", Mi, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, & IntegrandoReIntB0RehbergKlevanskymEnerA1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     }else{
//         // aux=1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskymkfinTfinRe", Mi, auxpolos,sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskymEnerB1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//         aux=aux+1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskymkfinTfinRe", Mi, auxpolos, sqrt(pow(Mi,2)+pow(Lambda,2)), &auxparams, &IntegrandoReIntB0RehbergKlevanskymEnerA1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
//     }   
//     return aux;
// }

double IntB0RehbergKlevanskypkfinTfinRe(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = T;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;

    vector<double> auxpolos={};
    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);

    vector<double> auxpolosp={};
    for (int i=0;i<auxpolos.size();i++){
        auxpolosp.push_back( sqrt(pow(auxpolos[i],2)-pow(Mi,2)));
    }

    aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskypkfinTfinRe", 0, Lambda, &auxparams, &IntegrandoReIntB0RehbergKlevanskyppB1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);

    if (auxpolos.size()==0){
        aux=aux+1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskypkfinTfinRe", 0, Lambda, &auxparams, &IntegrandoReIntB0RehbergKlevanskyppA1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }else{
        aux=aux+1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskypkfinTfinRe", 0, auxpolosp, Lambda, &auxparams, &IntegrandoReIntB0RehbergKlevanskyppA1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }   
    return aux;
}

double IntB0RehbergKlevanskymkfinTfinRe(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = T;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;
    
    vector<double> auxpolos={};
    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);
    
    vector<double> auxpolosp={};
    for (int i=0;i<auxpolos.size();i++){
        auxpolosp.push_back( sqrt(pow(auxpolos[i],2)-pow(Mi,2)));
    }

    aux=1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskymkfinTfinRe", 0, Lambda, &auxparams, &IntegrandoReIntB0RehbergKlevanskympB1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    if (auxpolos.size()==0){
        aux=aux+1./(2.*k)*integrate_QAGS_PRO("IntB0RehbergKlevanskymkfinTfinRe", 0, Lambda, &auxparams, & IntegrandoReIntB0RehbergKlevanskympA1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }else{
        aux=aux+1./(2.*k)*integrate_QAGP_PRO("IntB0RehbergKlevanskymkfinTfinRe", 0, auxpolosp, Lambda, &auxparams, &IntegrandoReIntB0RehbergKlevanskympA1D, int_Klev_workspace, f1_Klev_precision, f1_Klev_precision);
    }   
    return aux;
}

double IntB0RehbergKlevanskypkfinTfinIm(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = T;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;
    vector<double> auxpolos={};

    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);
    if (auxpolos.size()==0){
        aux=0;
    }else{
        if(auxpolos.size()==1){
            if(fabs(lambdax-(+k))>ZERO_Klev && fabs(lambdax-(-k))>ZERO_Klev){
                if(funEner1(lambdax, k, Mi, Mj)>=Mi && funEner1(lambdax, k, Mi, Mj)<=sqrt(pow(Mi,2)+pow(Lambda,2))){
                    aux=M_PI/k*doublesign(lambdax)*Intfermidist(T,funEner1(lambdax, k, Mi, Mj)-Cp,sqrt(pow(Mi,2)+pow(Lambda,2))-Cp);
                }else{
                    aux=M_PI/k*doublesign(lambdax)*Intfermidist(T,funEner2(lambdax, k, Mi, Mj)-Cp,sqrt(pow(Mi,2)+pow(Lambda,2))-Cp);
                }
            }else{
                if(fabs(Mi-Mj)>ZERO_Klev){
                    if(fabs(lambdax-(+k))<=ZERO_Klev){
                        aux=M_PI/k*doublesign(lambdax)*Intfermidist(T,+funEnerlambdaxkequal( k, Mi, Mj)-Cp,sqrt(pow(Mi,2)+pow(Lambda,2))-Cp);        
                    }else{
                        aux=M_PI/k*doublesign(lambdax)*Intfermidist(T, -funEnerlambdaxkequal( k, Mi, Mj)-Cp,sqrt(pow(Mi,2)+pow(Lambda,2))-Cp);
                    }
                }
            }
        }else{
            if(auxpolos.size()==2){
                aux=M_PI/k*doublesign(lambdax)*
                    Intfermidist(T,
                    min(funEner1(lambdax, k, Mi, Mj),funEner2(lambdax, k, Mi, Mj))-Cp,
                    max(funEner1(lambdax, k, Mi, Mj),funEner2(lambdax, k, Mi, Mj))-Cp);
            }
        }
    }
            
    return aux;
}

double IntB0RehbergKlevanskymkfinTfinIm(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;
    struct f1_loop_parameters auxparams;
    auxparams.temperature = T;
    auxparams.eff_chem_pot_1 = Cp;
    auxparams.cutoff = Lambda;
    auxparams.eff_mass_quark_1 = Mi;
    auxparams.eff_mass_quark_2 = Mj;
    auxparams.omega = lambdax;
    auxparams.external_momentum=k;
    vector<double> auxpolos={};

    auxpolos=polos(lambdax, Mi, Mj,  k, Lambda);

    if (auxpolos.size()==0){
        aux=0;
    }else{
        if(auxpolos.size()==1){
            if(fabs(lambdax-(+k))>ZERO_Klev && fabs(lambdax-(-k))>ZERO_Klev){
                if(funEner1(lambdax, k, Mi, Mj)>=Mi && funEner1(lambdax, k, Mi, Mj)<=sqrt(pow(Mi,2)+pow(Lambda,2))){
                    aux=-M_PI/k*doublesign(lambdax)*Intfermidist(T,-(funEner1(lambdax, k, Mi, Mj)+Cp),-(sqrt(pow(Mi,2)+pow(Lambda,2))+Cp));
                }else{
                    aux=-M_PI/k*doublesign(lambdax)*Intfermidist(T,-(funEner2(lambdax, k, Mi, Mj)+Cp),-(sqrt(pow(Mi,2)+pow(Lambda,2))+Cp));
                }
            }else{
                if(fabs(Mi-Mj)>ZERO_Klev){
                    if(fabs(lambdax-(+k))<=ZERO_Klev){
                        aux=-M_PI/k*doublesign(lambdax)*
                            Intfermidist(T,-(+funEnerlambdaxkequal( k, Mi, Mj)+Cp),-(sqrt(pow(Mi,2)+pow(Lambda,2))+Cp));
                    }else{
                        aux=-M_PI/k*doublesign(lambdax)*
                            Intfermidist(T,-(-funEnerlambdaxkequal( k, Mi, Mj)+Cp),-(sqrt(pow(Mi,2)+pow(Lambda,2))+Cp));
                    }
                }
            }
        }else{
            if(auxpolos.size()==2){
                aux=-M_PI/k*doublesign(lambdax)*
                    Intfermidist(T,
                    -(min(funEner1(lambdax, k, Mi, Mj),funEner2(lambdax, k, Mi, Mj))+Cp),
                    -(max(funEner1(lambdax, k, Mi, Mj),funEner2(lambdax, k, Mi, Mj))+Cp));
            }
        }
    }
            
    return aux;
}
//////////////////////
double IntB0RehbergKlevanskypRe(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;

    if(fabs(T)<=ZERO_Klev){
        if(fabs(k)<=ZERO_Klev)
        {
            aux=IntB0RehbergKlevanskypk0T0Re(lambdax,Mi,Mj,Lambda,Cp);
        }else{
            aux=IntB0RehbergKlevanskypkfinT0Re(lambdax,Mi,Mj,k,Lambda,Cp);
        }
    }else{
        if(fabs(k)<=ZERO_Klev){
            aux=IntB0RehbergKlevanskypk0TfinRe(lambdax,Mi,Mj,Lambda,T,Cp);
            }
        else{
            aux=IntB0RehbergKlevanskypkfinTfinRe(lambdax,Mi,Mj,k,Lambda,T, Cp);
        }
    }
    return aux;
}

double IntB0RehbergKlevanskymRe(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;

    if(fabs(T)<=ZERO_Klev){
        if(fabs(k)<=ZERO_Klev)
        {
            aux=IntB0RehbergKlevanskymk0T0Re(lambdax,Mi,Mj,Lambda,Cp);
        }else{
            aux=IntB0RehbergKlevanskymkfinT0Re(lambdax,Mi,Mj,k,Lambda,Cp);
        }
    }else{
        if(fabs(k)<=ZERO_Klev){
            aux=IntB0RehbergKlevanskymk0TfinRe(lambdax,Mi,Mj,Lambda,T,Cp);
            }
        else{
            aux=IntB0RehbergKlevanskymkfinTfinRe(lambdax,Mi,Mj,k,Lambda,T, Cp);   
        }
    }
    return aux;
}


double IntB0RehbergKlevanskypIm(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;
    
    if(fabs(T)<=ZERO_Klev){
        if(fabs(k)<=ZERO_Klev)
        {
            aux=IntB0RehbergKlevanskypk0T0Im(lambdax,Mi,Mj,Lambda,Cp);
        }else{
            aux=IntB0RehbergKlevanskypkfinT0Im(lambdax,Mi,Mj,k,Lambda,Cp);
        }
    }else{
        if(fabs(k)<=ZERO_Klev){
            aux=IntB0RehbergKlevanskypk0TfinIm(lambdax,Mi,Mj,Lambda,T,Cp);
            }
        else{
            aux=IntB0RehbergKlevanskypkfinTfinIm(lambdax, Mi, Mj, k, Lambda, T, Cp); 
        }
    }
    
    return aux;
}

double IntB0RehbergKlevanskymIm(double lambdax, double Mi, double Mj, double k, double Lambda, double T, double Cp){
    double aux=0;

    if(fabs(T)<=ZERO_Klev){
        if(fabs(k)<=ZERO_Klev)
        {
            aux=IntB0RehbergKlevanskymk0T0Im(lambdax,Mi,Mj,Lambda,Cp);
        }else{
            aux=IntB0RehbergKlevanskymkfinT0Im(lambdax,Mi,Mj,k,Lambda,Cp);
        }
    }else{
        if(fabs(k)<=ZERO_Klev){
            aux=IntB0RehbergKlevanskymk0TfinIm(lambdax,Mi,Mj,Lambda,T,Cp);
            }
        else{
            aux=IntB0RehbergKlevanskymkfinTfinIm(lambdax, Mi, Mj, k, Lambda, T, Cp); 
        }
    }
    return aux;
}


double IntB0RehbergKlevanskyk00k0Re(double Mi, double Mj, double Lambda, double T, double Cpi, double Cpj){
    double aux=real16Pi2f1Zero3Momentum3DCutoff(T, Cpi, Cpj, Lambda, Mi, Mj, 0.0, f1_Klev_precision);
    return aux;
}
double IntB0RehbergKlevanskyk00k0Im(double Mi, double Mj, double Lambda, double T, double Cpi, double Cpj){
    double aux=imag16Pi2f1Zero3Momentum3DCutoff(T, Cpi,Cpj, Lambda, Mi, Mj, 0.0);
    return aux;
}
double IntB0RehbergKlevanskyk0Re(double Mi, double Mj, double k0, double Lambda, double T, double Cpi, double Cpj){
    double aux=real16Pi2f1Zero3Momentum3DCutoff(T, Cpi, Cpj, Lambda, Mi, Mj, k0, f1_Klev_precision);
    return aux;
}
double IntB0RehbergKlevanskyk0Im(double Mi, double Mj, double k0, double Lambda, double T, double Cpi, double Cpj){
    double aux=imag16Pi2f1Zero3Momentum3DCutoff(T, Cpi,Cpj, Lambda, Mi, Mj, k0);
    return aux;
}

double IntB0RehbergKlevanskyRe(double Mi, double Mj, double k0,double k, double Lambda, double T, double Cpi, double Cpj){
    double aux=0;

    if( fabs(k0)<=ZERO_Klev && fabs(k)<=ZERO_Klev){
        aux=IntB0RehbergKlevanskyk00k0Re(Mi,Mj,Lambda,T,Cpi,Cpj);
    }else{
        if(fabs(k0)>ZERO_Klev && fabs(k)<=ZERO_Klev){
            aux=IntB0RehbergKlevanskyk0Re(Mi,Mj,k0,Lambda,T,Cpi,Cpj);
            }
        else{
            if((fabs(k0)<=ZERO_Klev&&fabs(k)>ZERO_Klev)||(fabs(k0)>ZERO_Klev&&fabs(k)>ZERO_Klev)){
                double lambdax= k0+Cpi-Cpj;
                aux=aux+IntB0RehbergKlevanskypRe(-lambdax,Mi,Mj,k,Lambda,T,Cpi);
                //cout<<"Re part1 "<<IntB0RehbergKlevanskypRe(-lambdax,Mi,Mj,k,Lambda,T,Cpi)<<"\n";                
                aux=aux-IntB0RehbergKlevanskymRe(+lambdax,Mi,Mj,k,Lambda,T,Cpi);
                //cout<<"Re part2 "<<IntB0RehbergKlevanskymRe(+lambdax,Mi,Mj,k,Lambda,T,Cpi)<<"\n";
                aux=aux+IntB0RehbergKlevanskypRe(+lambdax,Mj,Mi,k,Lambda,T,Cpj);
                //cout<<"Re part3 "<<IntB0RehbergKlevanskypRe(+lambdax,Mj,Mi,k,Lambda,T,Cpj)<<"\n";
                aux=aux-IntB0RehbergKlevanskymRe(-lambdax,Mj,Mi,k,Lambda,T,Cpj);
                //cout<<"Re part4 "<<IntB0RehbergKlevanskymRe(-lambdax,Mj,Mi,k,Lambda,T,Cpj)<<"\n";
                }
            }   
        }
    return aux;
}

double IntB0RehbergKlevanskyIm(double Mi, double Mj, double k0,double k, double Lambda, double T, double Cpi, double Cpj){
    double aux=0;

    if( fabs(k0)<=ZERO_Klev && fabs(k)<=ZERO_Klev){
        aux=IntB0RehbergKlevanskyk00k0Im(Mi,Mj,Lambda,T,Cpi,Cpj);
    }else{
        if(fabs(k0)>ZERO_Klev && fabs(k)<=ZERO_Klev){
            aux=IntB0RehbergKlevanskyk0Im(Mi,Mj,k0,Lambda,T,Cpi,Cpj);
            }
        else{
            if((fabs(k0)<=ZERO_Klev&&fabs(k)>ZERO_Klev)||(fabs(k0)>ZERO_Klev&&fabs(k)>ZERO_Klev)){
                double lambdax= k0+Cpi-Cpj;
                aux=aux+IntB0RehbergKlevanskypIm(-lambdax,Mi,Mj,k,Lambda,T,Cpi);
                //cout<<"Im part1 "<<IntB0RehbergKlevanskypIm(-lambdax,Mi,Mj,k,Lambda,T,Cpi)<<"\n";                
                aux=aux-IntB0RehbergKlevanskymIm(+lambdax,Mi,Mj,k,Lambda,T,Cpi);
                //cout<<"Im part2 "<<IntB0RehbergKlevanskymIm(+lambdax,Mi,Mj,k,Lambda,T,Cpi)<<"\n";
                aux=aux+IntB0RehbergKlevanskypIm(+lambdax,Mj,Mi,k,Lambda,T,Cpj);
                //cout<<"Im part3 "<<IntB0RehbergKlevanskypIm(+lambdax,Mj,Mi,k,Lambda,T,Cpj)<<"\n";
                aux=aux-IntB0RehbergKlevanskymIm(-lambdax,Mj,Mi,k,Lambda,T,Cpj);
                //cout<<"Im part4 "<<IntB0RehbergKlevanskymIm(-lambdax,Mj,Mi,k,Lambda,T,Cpj)<<"\n";
                }
            }   
        }
    return aux;
}

double Re16Pi2f1_Klev(double T, double mu1, double mu2, double Lambda, double M1, double M2, double w, double q)
{
    double Ref1_aux = IntB0RehbergKlevanskyRe(M1, M2, w, q, Lambda, T,  mu1, mu2);

    return Ref1_aux;
}

double Im16Pi2f1_Klev(double T, double mu1, double mu2, double Lambda, double M1, double M2, double w, double q)
{
    double Imf1_aux = IntB0RehbergKlevanskyIm(M1, M2, w, q, Lambda, T,  mu1, mu2);

    return Imf1_aux;
}

double Re16Pi2f1_Klev_Avg(double T, double mu1, double mu2, double Lambda, double M1, double M2, double w, double q)
{
    double Ref1_aux = 0.5*(IntB0RehbergKlevanskyRe(M1, M2, w, q, Lambda, T,  mu1, mu2)
                          +IntB0RehbergKlevanskyRe(M1, M2,-w, q, Lambda, T,  mu1, mu2));

    return Ref1_aux;
}

double Im16Pi2f1_Klev_Avg(double T, double mu1, double mu2, double Lambda, double M1, double M2, double w, double q)
{
    double Imf1_aux = 0.5*(IntB0RehbergKlevanskyIm(M1, M2, w, q, Lambda, T,  mu1, mu2)
                          +IntB0RehbergKlevanskyIm(M1, M2,-w, q, Lambda, T,  mu1, mu2));

    return Imf1_aux;
}

////////////////////////////////////////////////////////////

gsl_complex klevanskyB0Integral3DCutoffKlevanskyRecipe(
    NJL3DCutoffRegularizationScheme reguScheme, 
    double T, 
    double effCP1, 
    double effCP2, 
    double cutoff, 
    double M1, 
    double M2, 
    double w, 
    double k, 
    double integralPrecision
)
{
    if ( reguScheme==CUTOFF_ON_DIVERGENT_INTEGRALS_ONLY )
    {
        cout << "The function klevanskyB0Integral3DCutoff is not defined for the NJL3DCutoffRegularizationScheme:CUTOFF_ON_DIVERGENT_INTEGRALS_ONLY! Aborting!\n";
        abort();
    }

    double ReB0Klev = Re16Pi2f1_Klev_Avg(T, effCP1, effCP2, cutoff, M1, M2, w, k);
    double ImB0Klev = Im16Pi2f1_Klev_Avg(T, effCP1, effCP2, cutoff, M1, M2, w, k);

    gsl_complex B0Klev = gsl_complex_rect(ReB0Klev, ImB0Klev);

    return B0Klev;
}

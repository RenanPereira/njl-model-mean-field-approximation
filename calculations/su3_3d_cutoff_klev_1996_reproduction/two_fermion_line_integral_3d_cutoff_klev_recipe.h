#ifndef TWO_FERMION_LINE_INTEGRAL_3D_CUTOFF_KLEV_RECIPE_H
#define TWO_FERMION_LINE_INTEGRAL_3D_CUTOFF_KLEV_RECIPE_H

#include <gsl/gsl_complex_math.h>
#include "njl_model/njl_regularization_schemes.h"

struct f1_loop_parameters{
	double temperature;
	double eff_chem_pot_1;
	double eff_chem_pot_2;
	double cutoff;
	double eff_mass_quark_1;
	double eff_mass_quark_2;
	double omega;
	double external_momentum;
	double eta_variable;
};

double integral_inspector(std::string , void* , int);

double integrate_QAGS_PRO(std::string , double , double , void* , double (double , void*) , int , double , double );

double integrate_QAGP_PRO(std::string , double , std::vector<double> , double , void* , double (double , void*), int , double , double );

////////////////////////////////////////////////////////////

double Intfermidist(double , double , double );

double IntfermidistT0(double , double );

double doublesign(double );

double funEner1lambdax0(double , double , double );

double funEner1lambdaxnot0(double , double , double , double );

double funEner1(double , double , double , double );

double funEner2lambdax0(double , double , double );

double funEner2lambdaxnot0(double , double , double , double );

double funEner2(double , double , double , double );

double funEnerlambdaxkequal(double , double , double );

double IntegrandoReIntB0RehbergKlevanskypEner(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskymEner(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskypEnerA(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskypEnerB(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskymEnerA(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskymEnerB(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskypp(double , double , double , double , double , double , double , double);

double IntegrandoReIntB0RehbergKlevanskymp(double , double , double , double , double , double , double , double);

double IntegrandoReIntB0RehbergKlevanskyppA(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskyppB(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskympA(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskympB(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskypEner1D(double , void* );

double IntegrandoReIntB0RehbergKlevanskymEner1D(double , void* );

double IntegrandoReIntB0RehbergKlevanskypp1D(double , void* );

double IntegrandoReIntB0RehbergKlevanskymp1D(double , void* );

double IntB0RehbergKlevanskypkfinTfinRe(double , double , double , double , double , double , double );

double IntB0RehbergKlevanskymkfinTfinRe(double , double , double , double , double , double , double );

double IntB0RehbergKlevanskypkfinTfinReAlt(double , double , double , double , double , double , double );

double IntB0RehbergKlevanskymkfinTfinReAlt(double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskypEnerA1D(double , void*);

double IntegrandoReIntB0RehbergKlevanskypEnerB1D(double , void*);

double IntegrandoReIntB0RehbergKlevanskymEnerA1D(double , void*);

double IntegrandoReIntB0RehbergKlevanskymEnerB1D(double , void*);

double IntegrandoReIntB0RehbergKlevanskyppA1D(double , void*);

double IntegrandoReIntB0RehbergKlevanskyppB1D(double , void*);

double IntegrandoReIntB0RehbergKlevanskympA1D(double , void*);

double IntegrandoReIntB0RehbergKlevanskympB1D(double , void*);

double IntegrandoReIntB0RehbergKlevanskypT0Ener(double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskymT0Ener(double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskypT0Enerlambdax0(double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskymT0Enerlambdax0(double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskypT0Ener2(double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskymT0Ener2(double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskypT0Ener21D(double , void* );

double IntegrandoReIntB0RehbergKlevanskymT0Ener21D(double , void* );

std::vector<double> polos(double , double , double , double , double );

double IntB0RehbergKlevanskypk0T0Re(double , double , double , double , double );

double IntB0RehbergKlevanskypk0T0Im(double , double , double , double , double );

double IntB0RehbergKlevanskymk0T0Re(double , double , double , double , double );

double IntB0RehbergKlevanskymk0T0Im(double , double , double , double , double );

double IntB0RehbergKlevanskypk0TfinRe(double , double , double , double , double , double );

double IntB0RehbergKlevanskypk0TfinIm(double , double , double , double , double , double );

double IntB0RehbergKlevanskymk0TfinRe(double , double , double , double , double , double );

double IntB0RehbergKlevanskymk0TfinIm(double , double , double , double , double , double );

double IntB0RehbergKlevanskypkfinT0Re(double , double , double , double , double , double );

double IntB0RehbergKlevanskymkfinT0Re(double , double , double , double , double , double );

double IntB0RehbergKlevanskypkfinT0Im(double , double , double , double , double , double );

double IntB0RehbergKlevanskymkfinT0Im(double , double , double , double , double , double );

double IntB0RehbergKlevanskypkfinTfinIm(double , double , double , double , double , double , double );

double IntB0RehbergKlevanskymkfinTfinIm(double , double , double , double , double , double , double );

double IntB0RehbergKlevanskypRe(double , double , double , double , double , double , double);

double IntB0RehbergKlevanskymRe(double , double , double , double , double , double , double);

double IntB0RehbergKlevanskypIm(double , double , double , double , double , double , double);

double IntB0RehbergKlevanskymIm(double , double , double , double , double , double , double);

double IntB0RehbergKlevanskyRe(double , double , double ,double , double , double , double , double);

double IntB0RehbergKlevanskyIm(double , double , double ,double , double , double , double , double);

double IntB0RehbergKlevanskyk00k0Re(double , double , double , double , double , double );

double IntB0RehbergKlevanskyk00k0Im(double , double , double , double , double , double );

double IntB0RehbergKlevanskyk0Re(double , double , double , double , double , double , double );

double IntB0RehbergKlevanskyk0Im(double , double , double , double , double , double , double );

double Re16Pi2f1_Klev(double , double , double , double , double , double , double , double );

double Im16Pi2f1_Klev(double , double , double , double , double , double , double , double );

double Re16Pi2f1_Klev_Avg(double , double , double , double , double , double , double , double );

double Im16Pi2f1_Klev_Avg(double , double , double , double , double , double , double , double );

gsl_complex klevanskyB0Integral3DCutoffKlevanskyRecipe(
    NJL3DCutoffRegularizationScheme  , 
    double , 
    double , 
    double , 
    double , 
    double , 
    double , 
    double , 
    double , 
    double 
);

#endif
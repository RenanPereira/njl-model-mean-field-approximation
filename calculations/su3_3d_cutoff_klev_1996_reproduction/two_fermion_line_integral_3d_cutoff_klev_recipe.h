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

double IntegrandoReIntB0RehbergKlevanskypEner(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskymEner(double , double , double , double , double , double , double , double );

double IntegrandoReIntB0RehbergKlevanskypp(double , double , double , double , double , double , double , double);

double IntegrandoReIntB0RehbergKlevanskymp(double , double , double , double , double , double , double , double);

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

double IntB0RehbergKlevanskypRe(double , double , double , double , double , double , double);

double IntB0RehbergKlevanskymRe(double , double , double , double , double , double , double);

double IntB0RehbergKlevanskypIm(double , double , double , double , double , double , double);

double IntB0RehbergKlevanskymIm(double , double , double , double , double , double , double);

double IntB0RehbergKlevanskyRe(double , double , double ,double , double , double , double , double);

double IntB0RehbergKlevanskyIm(double , double , double ,double , double , double , double , double);

double IntB0RehbergKlevanskyk00k0Re(double , double , double , double , double , double );

double IntB0RehbergKlevanskyk00k0Im(double , double , double , double , double , double );

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
#ifndef N_FERMION_LINE_INTEGRALS_CALCULATOR_H
#define N_FERMION_LINE_INTEGRALS_CALCULATOR_H

#include <string>
#include "njl_model/njl_regularization_schemes.h"

std::string trim0ToDot0(const double );

void evaluateKlevanskyB0Integral3DCutoffVsZeroMomentumToFile(
    const int , 
    const double , 
    const double ,
    const NJL3DCutoffRegularizationScheme , 
    const double ,
    const double ,
    const double ,
    const double ,
    const double ,
    const double ,
    const double ,
    const double );

void evaluateKlevanskyB0Integral3DCutoffVsThreeMomentumToFile(
    const int , 
    const double , 
    const double ,
    const NJL3DCutoffRegularizationScheme , 
    const double ,
    const double ,
    const double ,
    const double ,
    const double ,
    const double ,
    const double ,
    const double );

#endif
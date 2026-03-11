#ifndef N_FERMION_LINE_INTEGRALS_CALCULATOR_H
#define N_FERMION_LINE_INTEGRALS_CALCULATOR_H

#include "njl_model/njl_regularization_schemes.h"

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
    const double 
);

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
    const double 
);

#endif
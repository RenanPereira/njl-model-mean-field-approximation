#ifndef ONE_FERMION_LINE_INTEGRAL_3D_CUTOFF_H
#define ONE_FERMION_LINE_INTEGRAL_3D_CUTOFF_H

#include <iostream>
#include <gsl/gsl_complex.h>
#include "gsl_wrapper/Integration1DimGSL.h"
#include "njl_model/njl_regularization_schemes.h"


class OneFermionLine3DCutoffIntegrand : public GeneralIntegrandParameters
{
public:
    std::string integralID;
    double temperature = 0.0/0.0;
    double effectiveChemicalPotential = 0.0/0.0;
    double effectiveMass = 0.0/0.0;

public:
    OneFermionLine3DCutoffIntegrand(std::string integralIDAux, double temperatureAux, double effectiveChemicalPotentialAux, double effectiveMassAux)
    {   
        integralID = integralIDAux;
        temperature = temperatureAux;
        effectiveChemicalPotential = effectiveChemicalPotentialAux;
        effectiveMass = effectiveMassAux;
    };

    OneFermionLine3DCutoffIntegrand(void* auxiliar)
    {   
        integralID = ((class OneFermionLine3DCutoffIntegrand *)(auxiliar))->integralID;
        temperature = ((class OneFermionLine3DCutoffIntegrand *)(auxiliar))->temperature;
        effectiveChemicalPotential = ((class OneFermionLine3DCutoffIntegrand *)(auxiliar))->effectiveChemicalPotential;
        effectiveMass = ((class OneFermionLine3DCutoffIntegrand *)(auxiliar))->effectiveMass;
    };

    std::string getIntegralID(){ return integralID; }
    double getTemperature(){ return temperature; };
    double getEffectiveChemicalPotential(){ return effectiveChemicalPotential; };
    double getEffectiveMass(){ return effectiveMass; };

    void printIntegrandVariables() override
    {   
        std::cout << integralID << "\n";
        std::cout << "T = " << temperature << "\n";
        std::cout << "effChemPot = " << effectiveChemicalPotential << "\n";
        std::cout << "effMass = " << effectiveMass << "\n";
        std::cout << "(T, effChemPot, effMass)\n";
        std::cout << "(" << temperature << ", " << effectiveChemicalPotential << ", " << effectiveMass << ")" << "\n";
    }
};

double f0Primitive(double , double );

double f0ConvergentIntegrand(double , void *);

double f03DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double );

double realKlevanskyAIntegral3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double );

gsl_complex klevanskyAIntegral3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double );

double sigmaNJL3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double , double );

double gEEta(double , double , double , double , double , double , double , double );

double realKlevanskyAPair3DCutoffIntegrand(double , void *);

double realKlevanskyAPair3DCutoffM2LargerM1(double , double , double , double , double , double , double );

double realKlevanskyAPair3DCutoffM2EqualM1(double , double , double , double , double , double , double );

double realKlevanskyAPair3DCutoff(double , double , double , double , double , double , double );

double gEpsilonEta(double , double , double , double , double , double , double , double );

double realKlevanskyAScat3DCutoffIntegrand(double , void *);

double realKlevanskyAScat3DCutoffM2LargerM1(double , double , double , double , double , double , double );

double realKlevanskyAScat3DCutoffM2EqualM1(double , double , double , double , double , double , double );

double realKlevanskyAScat3DCutoff(double , double , double , double , double , double , double );

double realKlevanskyAIntegral3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double , double );

gsl_complex klevanskyAIntegral3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double , double );


#endif
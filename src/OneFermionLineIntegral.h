#ifndef ONEFERMIONLINEINTEGRAL_H
#define ONEFERMIONLINEINTEGRAL_H

#include <iostream>
#include <gsl/gsl_complex.h>
#include "gsl_wrapper/Integration1DimGSL.h"
#include "generalPhysicsAndMath.h"


class OneFermionLine3DCutoffIntegrand : public GeneralIntegrandParameters
{
public:
    string integralID;
    double temperature = 0.0/0.0;
    double effectiveChemicalPotential = 0.0/0.0;
    double effectiveMass = 0.0/0.0;

public:
    OneFermionLine3DCutoffIntegrand(string integralIDAux, double temperatureAux, double effectiveChemicalPotentialAux, double effectiveMassAux)
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

    string getIntegralID(){ return integralID; }
    double getTemperature(){ return temperature; };
    double getEffectiveChemicalPotential(){ return effectiveChemicalPotential; };
    double getEffectiveMass(){ return effectiveMass; };

    void printIntegrandVariables() override
    {   
        cout << integralID << "\n";
        cout << "T = " << temperature << "\n";
        cout << "effChemPot = " << effectiveChemicalPotential << "\n";
        cout << "effMass = " << effectiveMass << "\n";
        cout << "(T, effChemPot, effMass)\n";
        cout << "(" << temperature << ", " << effectiveChemicalPotential << ", " << effectiveMass << ")" << "\n";
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
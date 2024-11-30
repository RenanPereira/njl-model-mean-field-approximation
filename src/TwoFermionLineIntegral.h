#ifndef TWOFERMIONLINEINTEGRAL_H
#define TWOFERMIONLINEINTEGRAL_H

#include <iostream>
#include <gsl/gsl_complex.h>
#include "gsl_wrapper/Integration1DimGSL.h"
#include "generalPhysicsAndMath.h"


class TwoFermionLine3DCutoffIntegrand : public GeneralIntegrandParameters
{
private:
    string integralID = "notDefined";
    double threeMomentumCutoff = 0.0/0.0;
    double temperature = 0.0/0.0;
    double effectiveChemicalPotential1 = 0.0/0.0;
    double effectiveChemicalPotential2 = 0.0/0.0;
    double effectiveMass1 = 0.0/0.0;
    double effectiveMass2 = 0.0/0.0;
    double zeroMomentum = 0.0/0.0;
    double threeMomentum = 0.0/0.0;
    double etaVariable = 0.0/0.0;

public:
    TwoFermionLine3DCutoffIntegrand(string integralIDAux, 
                                    double threeMomentumCutoffAux, 
                                    double effectiveMass1Aux, double effectiveMass2Aux, 
                                    double threeMomentumAux)
    {   
        integralID = integralIDAux;
        threeMomentumCutoff = threeMomentumCutoffAux;
        effectiveMass1 = effectiveMass1Aux;
        effectiveMass2 = effectiveMass2Aux;
        threeMomentum = threeMomentumAux;
    };

    TwoFermionLine3DCutoffIntegrand(string integralIDAux, 
                                    double temperatureAux, 
                                    double effectiveChemicalPotential1Aux, double effectiveChemicalPotential2Aux,
                                    double threeMomentumCutoffAux,
                                    double effectiveMass1Aux, double effectiveMass2Aux,
                                    double zeroMomentumAux, double threeMomentumAux,
                                    double etaVariableAux)
    {   
        integralID = integralIDAux;
        temperature = temperatureAux;
        effectiveChemicalPotential1 = effectiveChemicalPotential1Aux;
        effectiveChemicalPotential2 = effectiveChemicalPotential2Aux;
        threeMomentumCutoff = threeMomentumCutoffAux;
        effectiveMass1 = effectiveMass1Aux;
        effectiveMass2 = effectiveMass2Aux;
        zeroMomentum = zeroMomentumAux;
        threeMomentum = threeMomentumAux;
        etaVariable = etaVariableAux;
    };

    TwoFermionLine3DCutoffIntegrand(double temperatureAux, 
                                    double effectiveChemicalPotential1Aux, double effectiveChemicalPotential2Aux,
                                    double threeMomentumCutoffAux,
                                    double effectiveMass1Aux, double effectiveMass2Aux,
                                    double zeroMomentumAux, double threeMomentumAux)
    {   
        temperature = temperatureAux;
        effectiveChemicalPotential1 = effectiveChemicalPotential1Aux;
        effectiveChemicalPotential2 = effectiveChemicalPotential2Aux;
        threeMomentumCutoff = threeMomentumCutoffAux;
        effectiveMass1 = effectiveMass1Aux;
        effectiveMass2 = effectiveMass2Aux;
        zeroMomentum = zeroMomentumAux;
        threeMomentum = threeMomentumAux;
    };

    TwoFermionLine3DCutoffIntegrand(void* auxiliar)
    {   
        integralID = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->integralID;
        threeMomentumCutoff = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->threeMomentumCutoff;
        temperature = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->temperature;
        effectiveChemicalPotential1 = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->effectiveChemicalPotential1;
        effectiveChemicalPotential2 = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->effectiveChemicalPotential2;
        effectiveMass1 = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->effectiveMass1;
        effectiveMass2 = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->effectiveMass2;
        zeroMomentum = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->zeroMomentum;
        threeMomentum = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->threeMomentum;
        etaVariable = ((class TwoFermionLine3DCutoffIntegrand *)(auxiliar))->etaVariable;
    };

    string getIntegralID(){ return integralID; }
    double getThreeMomentumCutoff(){ return threeMomentumCutoff; };
    double getTemperature(){ return temperature; };
    double getEffectiveChemicalPotential1(){ return effectiveChemicalPotential1; };
    double getEffectiveChemicalPotential2(){ return effectiveChemicalPotential2; };
    double getEffectiveMass1(){ return effectiveMass1; };
    double getEffectiveMass2(){ return effectiveMass2; };
    double getZeroMomentum(){ return zeroMomentum; };
    double getThreeMomentum(){ return threeMomentum; };
    double getEtaVariable(){ return etaVariable; };

    void setIntegralID(string integralIDAux){ integralID = integralIDAux; }
    void setEtaVariable(double etaVariableAux){ etaVariable = etaVariableAux; };

    void printIntegrandVariables() override
    {   
        printf("integralID = %s\n", integralID.c_str());
        printf("cutoff3D = %.20f\n", threeMomentumCutoff);
        printf("T = %.20f\n", temperature);
        printf("effChemPot1 = %.20f\n", effectiveChemicalPotential1);
        printf("effChemPot2 = %.20f\n", effectiveChemicalPotential2);
        printf("effMass1 = %.20f\n", effectiveMass1);
        printf("effMass2 = %.20f\n", effectiveMass2);
        printf("zeroMom = %.20f\n", zeroMomentum);
        printf("threeMom = %.20f\n", threeMomentum);
        printf("eta = %.20f\n", etaVariable);
        /*
        int defaultPrecision = int(cout.precision());
        cout << setprecision(15); //must include #include <iomanip>
        cout << "integralID = " << integralID << "\n";
        cout << "cutoff3D = " << threeMomentumCutoff << "\n";
        cout << "T = " << temperature << "\n";
        cout << "effChemPot1 = " << effectiveChemicalPotential1 << "\n";
        cout << "effChemPot2 = " << effectiveChemicalPotential2 << "\n";
        cout << "effMass1 = " << effectiveMass1 << "\n";
        cout << "effMass2 = " << effectiveMass2 << "\n";
        cout << "zeroMom = " << zeroMomentum << "\n";
        cout << "threeMom = " << threeMomentum << "\n";
        cout << "eta = "  << etaVariable << "\n";
        cout << setprecision(defaultPrecision);
        cout << "(cutoff3D, T, effChemPot1, effChemPot2, effMass1, effMass2, extZeroMom, extThreeMom, eta)\n";
        cout << "(" << threeMomentumCutoff << ", " << temperature << ", " 
                    << effectiveChemicalPotential1 << ", " << effectiveChemicalPotential2 << ", "
                    << effectiveMass1 << ", " << effectiveMass2 << ", "
                    << zeroMomentum << ", " << threeMomentum << ", "
                    << etaVariable << ")" << "\n";
        */
    }
};


double imag16Pi2f1PairSign(double , bool );

double imag16Pi2f1ScatSign(double , double , bool );

double dEdepsilon_epsilonMinPlus(double , double , double , double );

double dEdepsilon_epsilonMinMinus(double , double , double , double );

double dEdepsilon_epsilonLambdaM1Plus(double , double , double );

double dEdepsilon_epsilonLambdaM1Minus(double , double , double );

double dEdepsilon_epsilonLambdaM2Plus(double , double , double );

double dEdepsilon_epsilonLambdaM2Minus(double , double , double );

double dEdepsilon_gamma1E(double , double , double , double );

double dEdepsilon_gamma1epsilon(double , double , double , double );

double dEdepsilon_gamma2E(double , double , double , double );

double dEdepsilon_gamma2epsilon(double , double , double , double );

double dEdepsilon_alpha1E(double , double , double );

double dEdepsilon_alpha1epsilon(double , double , double );

double dEdepsilon_LambdaSwitchE(double , double , double );

double dEdepsilon_LambdaSwitchepsilon(double , double , double );

double dEdepsilon_epsilonMax(double , double , double , double , double );

double dEdepsilon_epsilonMin(double , double , double , double , double );

double dEdepsilon_EMinM1LargerM2(double , double , double , double );

double dEdepsilon_EMinM2LargerM1(double , double , double , double );

double dEdepsilon_EMin(double , double , double , double );

double dEdepsilon_EMax(double , double , double , double );

double dEdepsilon_AreaIntegrand(double , void *);

double dEdepsilon_Area(double , double , double , double , double );

double depsilondE_EMinPlus(double , double , double , double );

double depsilondE_EMinMinus(double , double , double , double );

double depsilondE_ELambdaM1Plus(double , double , double );

double depsilondE_ELambdaM1Minus(double , double , double );

double depsilondE_ELambdaM2Plus(double , double , double );

double depsilondE_ELambdaM2Minus(double , double , double );

double depsilondE_LambdaSwitchepsilon(double , double , double );

double depsilondE_ELambdaPlus(double , double , double , double );

double depsilondE_gamma1epsilon(double , double , double , double );

double depsilondE_gamma2epsilon(double , double , double , double );

double depsilondE_gamma1E(double , double , double , double );

double depsilondE_gamma2E(double , double , double , double );

double depsilondE_alpha1epsilon(double , double , double );

double depsilondE_alpha2epsilon(double , double , double );

double depsilondE_alpha1E(double , double , double );

double depsilondE_alpha2E(double , double , double );

double depsilondE_EMaxM1LargerM2(double , double , double , double , double );

double depsilondE_EMaxM2LargerM1(double , double , double , double , double );

double depsilondE_EMax(double , double , double , double , double );

double depsilondE_EMin(double , double , double , double );

double depsilondE_epsilonMin(double , double , double , double );

double depsilondE_epsilonMax(double , double , double , double );

double depsilondE_AreaIntegrand(double , void *);

double depsilondE_Area(double , double , double , double , double );

double gPlusEta(double , double , double , double , double , double , double , double , double );

double real16Pi2f1Pair3DCutoffNumerator(double , void *);

double real16Pi2f1Pair3DCutoff(double , double , double , double , double , double , double , double , double );

double imag16Pi2f1Pair3DCutoff(double , double , double , double , double , double , double , double );

double gMinusEta(double , double , double , double , double , double , double , double , double );

double real16Pi2f1Scat3DCutoffNumerator(double , void *);

double real16Pi2f1Scat3DCutoff(double , double , double , double , double , double , double , double , double );

double imag16Pi2f1Scat3DCutoff(double , double , double , double , double , double , double , double );

double real16Pi2f1Finite3Momentum3DCutoff(double , double , double , double , double , double , double , double , double );

double imag16Pi2f1Finite3Momentum3DCutoff(double , double , double , double , double , double , double , double );

double E0(double , double );

double ELambda(double , double , double );

double pFunctionE(double , double , double );

double E1FunctionE(double , double , double );

double E2FunctionE(double , double , double );

double gPlusEtaZero3Momentum(double , double , double , double , double , double , double );

double real16Pi2f1PairZero3Momentum3DCutoffNumerator(double , void *);

double real16Pi2f1PairZero3Momentum3DCutoff(double , double , double , double , double , double , double , double );

double imag16Pi2f1PairZero3Momentum3DCutoff(double , double , double , double , double , double , double );

double epsilon0(double , double );

double epsilonLambda(double , double , double );

double pFunctionEpsilon(double , double , double );

double E1FunctionEpsilon(double , double , double );

double E2FunctionEpsilon(double , double , double );

double gMinusEtaZero3Momentum(double , double , double , double , double , double , double );

double gMinusEtaZero3Momentum(double , double , double , double , double , double );

double real16Pi2f1ScatZero3Momentum3DCutoffDifferentMassesNumerator(double , void *);

double real16Pi2f1ScatZero3Momentum3DCutoffDifferentMasses(double , double , double , double , double , double , double , double );

double real16Pi2f1ScatZero3Momentum3DCutoffEqualMassesIntegrand(double , void *);

double real16Pi2f1ScatZero3Momentum3DCutoffEqualMasses(double , double , double , double , double , double , double , double );

double imag16Pi2f1ScatZero3Momentum3DCutoffDifferentMasses(double , double , double , double , double , double , double );

double real16Pi2f1Zero3Momentum3DCutoff(double , double , double , double , double , double , double , double );

double imag16Pi2f1Zero3Momentum3DCutoff(double , double , double , double , double , double , double );

double real16Pi2f13DCutoff(double , double , double , double , double , double , double , double , double );

double imag16Pi2f13DCutoff(double , double , double , double , double , double , double , double );

gsl_complex klevanskyB0Integral3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double , double , double , double , double );


#endif
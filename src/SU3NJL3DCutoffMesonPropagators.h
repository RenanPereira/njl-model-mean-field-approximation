#ifndef SU3NJL3DCUTOFFMESONPROPAGATORS_H
#define SU3NJL3DCUTOFFMESONPROPAGATORS_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "ComplexSquareMatrixGSL.h"
#include "rootSolverGSL.h"
#include "SU3NJL3DCutoff.h"


gsl_complex pseudoscalarPolarizationOperator3DCutoff(NJL3DCutoffRegularizationScheme , double , double , 
                                                     double , double , double , double , double , 
                                                     double , double , double , double );

gsl_complex pseudoscalarPolarizationOperator3DCutoff(NJL3DCutoffRegularizationScheme , double , double , 
                                                     double , double , double , double , double , 
                                                     double , double , double );


gsl_complex scalarPolarizationOperator3DCutoff(NJL3DCutoffRegularizationScheme , double , double , 
                                               double , double , double , double , double , 
                                               double , double , double , double );

gsl_complex scalarPolarizationOperator3DCutoff(NJL3DCutoffRegularizationScheme , double , double , 
                                               double , double , double , double , double , 
                                               double , double , double );


gsl_complex pionPlusPropagator(SU3NJL3DCutoffParameters , double , 
                               double , double , double , 
                               double , double , double , 
                               double , double , double , double );

gsl_complex pionPlusPropagator(SU3NJL3DCutoffParameters , double , 
                               double , double , double , 
                               double , double , double , 
                               double , double , double );

gsl_complex pionMinusPropagator(SU3NJL3DCutoffParameters , double , 
                                double , double , double , 
                                double , double , double , 
                                double , double , double , double );

gsl_complex pionMinusPropagator(SU3NJL3DCutoffParameters , double , 
                                double , double , double , 
                                double , double , double , 
                                double , double , double );

gsl_complex kaonPlusPropagator(SU3NJL3DCutoffParameters , double , 
                               double , double , double , 
                               double , double , double , 
                               double , double , double , double );

gsl_complex kaonPlusPropagator(SU3NJL3DCutoffParameters , double , 
                               double , double , double , 
                               double , double , double , 
                               double , double , double );

gsl_complex kaonMinusPropagator(SU3NJL3DCutoffParameters , double , 
                                double , double , double , 
                                double , double , double , 
                                double , double , double , double );

gsl_complex kaonMinusPropagator(SU3NJL3DCutoffParameters , double , 
                                double , double , double , 
                                double , double , double , 
                                double , double , double );

gsl_complex neutralKaonPropagator(SU3NJL3DCutoffParameters , double , 
                                  double , double , double , 
                                  double , double , double , 
                                  double , double , double , double );

gsl_complex neutralKaonPropagator(SU3NJL3DCutoffParameters , double , 
                                  double , double , double , 
                                  double , double , double , 
                                  double , double , double );

gsl_complex antiNeutralKaonPropagator(SU3NJL3DCutoffParameters , double , 
                                      double , double , double , 
                                      double , double , double , 
                                      double , double , double , double );

gsl_complex antiNeutralKaonPropagator(SU3NJL3DCutoffParameters , double , 
                                      double , double , double , 
                                      double , double , double , 
                                      double , double , double );

ComplexSquareMatrixGSL neutral038PseudoscalarsPropagator(SU3NJL3DCutoffParameters , double , 
                                                         double , double , double , 
                                                         double , double , double , 
                                                         double , double , double , double );

ComplexSquareMatrixGSL neutral038PseudoscalarsPropagator(SU3NJL3DCutoffParameters , double , 
                                                         double , double , double , 
                                                         double , double , double , 
                                                         double , double , double );

gsl_complex sigmaPionPlusPropagator(SU3NJL3DCutoffParameters , double , 
                                    double , double , double , 
                                    double , double , double , 
                                    double , double , double , double );

gsl_complex sigmaPionPlusPropagator(SU3NJL3DCutoffParameters , double , 
                                    double , double , double , 
                                    double , double , double , 
                                    double , double , double );

gsl_complex sigmaPionMinusPropagator(SU3NJL3DCutoffParameters , double , 
                                     double , double , double , 
                                     double , double , double , 
                                     double , double , double , double );

gsl_complex sigmaPionMinusPropagator(SU3NJL3DCutoffParameters , double , 
                                     double , double , double , 
                                     double , double , double , 
                                     double , double , double );

gsl_complex sigmaKaonPlusPropagator(SU3NJL3DCutoffParameters , double , 
                                    double , double , double , 
                                    double , double , double , 
                                    double , double , double , double );

gsl_complex sigmaKaonPlusPropagator(SU3NJL3DCutoffParameters , double , 
                                    double , double , double , 
                                    double , double , double , 
                                    double , double , double );

gsl_complex sigmaKaonMinusPropagator(SU3NJL3DCutoffParameters , double , 
                                     double , double , double , 
                                     double , double , double , 
                                     double , double , double , double );

gsl_complex sigmaKaonMinusPropagator(SU3NJL3DCutoffParameters , double , 
                                     double , double , double , 
                                     double , double , double , 
                                     double , double , double );

gsl_complex neutralSigmaKaonPropagator(SU3NJL3DCutoffParameters , double , 
                                       double , double , double , 
                                       double , double , double , 
                                       double , double , double , double );

gsl_complex neutralSigmaKaonPropagator(SU3NJL3DCutoffParameters , double , 
                                       double , double , double , 
                                       double , double , double , 
                                       double , double , double );

gsl_complex antiNeutralSigmaKaonPropagator(SU3NJL3DCutoffParameters , double , 
                                           double , double , double , 
                                           double , double , double , 
                                           double , double , double , double );

gsl_complex antiNeutralSigmaKaonPropagator(SU3NJL3DCutoffParameters , double , 
                                           double , double , double , 
                                           double , double , double , 
                                           double , double , double );

ComplexSquareMatrixGSL neutral038ScalarsPropagator(SU3NJL3DCutoffParameters , double , 
                                                   double , double , double , 
                                                   double , double , double , 
                                                   double , double , double , double );

ComplexSquareMatrixGSL neutral038ScalarsPropagator(SU3NJL3DCutoffParameters , double , 
                                                   double , double , double , 
                                                   double , double , double , 
                                                   double , double , double );


enum mesonState { pionPlus, pionMinus, kaonPlus, kaonMinus, neutralKaon, antiNeutralKaon, diagonalPseudoscalars,
                  sigmaPionPlus, sigmaPionMinus, sigmaKaonPlus, sigmaKaonMinus, neutralSigmaKaon, antiNeutralSigmaKaon, diagonalScalars };


gsl_complex nonDiagonalMesonPropagator(SU3NJL3DCutoffParameters , double , 
                                       double , double , double , 
                                       double , double , double , 
                                       double , double , double , double ,
                                       mesonState );


gsl_complex nonDiagonalMesonPropagator(SU3NJL3DCutoffParameters , double , 
                                       double , double , double , 
                                       double , double , double , 
                                       double , double , double ,
                                       mesonState );


ComplexSquareMatrixGSL diagonalMesonPropagator(SU3NJL3DCutoffParameters , double , 
                                               double , double , double , 
                                               double , double , double , 
                                               double , double , double , double ,
                                               mesonState );


ComplexSquareMatrixGSL diagonalMesonPropagator(SU3NJL3DCutoffParameters , double , 
                                               double , double , double , 
                                               double , double , double , 
                                               double , double , double ,
                                               mesonState );


int SU3NJL3DCutoffMesonMassEquations(const gsl_vector *, void *, gsl_vector *);


class SU3NJL3DCutoffMeson
{
private:
    SU3NJL3DCutoffParameters parametersNJL;
    double temperature;
    double upQuarkEffectiveChemicalPotential;
    double downQuarkEffectiveChemicalPotential;
    double strangeQuarkEffectiveChemicalPotential;
    double upQuarkEffectiveMass; 
    double downQuarkEffectiveMass;
    double strangeQuarkEffectiveMass;
    double integralPrecision;
    mesonState meson;

    double mesonMass = 0;
    double mesonWidth = 0;

public:
    SU3NJL3DCutoffMeson(SU3NJL3DCutoffParameters parametersNJLAux, double TAux, 
                        double effChemPotUAux, double effChemPotDAux, double effChemPotSAux, 
                        double effMassUAux, double effMassDAux, double effMassSAux, 
                        double integralPrecisionAux, mesonState mesonAux)
    {
        parametersNJL = parametersNJLAux;
        temperature = TAux;
        upQuarkEffectiveChemicalPotential = effChemPotUAux;
        downQuarkEffectiveChemicalPotential = effChemPotDAux;
        strangeQuarkEffectiveChemicalPotential = effChemPotSAux;
        upQuarkEffectiveMass = effMassUAux; 
        downQuarkEffectiveMass = effMassDAux;
        strangeQuarkEffectiveMass = effMassSAux;
        integralPrecision = integralPrecisionAux;
        meson = mesonAux;
    };

    SU3NJL3DCutoffMeson(void* auxiliar)
    {   
        parametersNJL = ((class SU3NJL3DCutoffMeson *)(auxiliar))->parametersNJL;
        temperature = ((class SU3NJL3DCutoffMeson *)(auxiliar))->temperature;
        upQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffMeson *)(auxiliar))->upQuarkEffectiveChemicalPotential;
        downQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffMeson *)(auxiliar))->downQuarkEffectiveChemicalPotential;
        strangeQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffMeson *)(auxiliar))->strangeQuarkEffectiveChemicalPotential;
        upQuarkEffectiveMass = ((class SU3NJL3DCutoffMeson *)(auxiliar))->upQuarkEffectiveMass;
        downQuarkEffectiveMass = ((class SU3NJL3DCutoffMeson *)(auxiliar))->downQuarkEffectiveMass;
        strangeQuarkEffectiveMass = ((class SU3NJL3DCutoffMeson *)(auxiliar))->strangeQuarkEffectiveMass;
        integralPrecision = ((class SU3NJL3DCutoffMeson *)(auxiliar))->integralPrecision;
        meson = ((class SU3NJL3DCutoffMeson *)(auxiliar))->meson;
    };

    SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };
    double getTemperature(){ return temperature; };
    double getUpQuarkEffectiveChemicalPotential(){ return upQuarkEffectiveChemicalPotential; };
    double getDownQuarkEffectiveChemicalPotential(){ return downQuarkEffectiveChemicalPotential; };
    double getStrangeQuarkEffectiveChemicalPotential(){ return strangeQuarkEffectiveChemicalPotential; };
    double getUpQuarkEffectiveMass(){ return upQuarkEffectiveMass; };
    double getDownQuarkEffectiveMass(){ return downQuarkEffectiveMass; };
    double getStrangeQuarkEffectiveMass(){ return strangeQuarkEffectiveMass; };
    double getMesonPropagatorIntegralPrecision(){ return integralPrecision; };
    mesonState getMesonState(){ return meson; }

    double getMesonMass(){ return mesonMass; }
    double getMesonWidth(){ return mesonWidth; }

    gsl_complex calculateNonDiagonalPropagator(double zeroMomentum, double threeMomentum, double gamma)
    {   
        gsl_complex propagator;
        propagator = nonDiagonalMesonPropagator(parametersNJL, temperature, 
                                                upQuarkEffectiveChemicalPotential, downQuarkEffectiveChemicalPotential, strangeQuarkEffectiveChemicalPotential, 
                                                upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                                                zeroMomentum, threeMomentum, gamma, integralPrecision,
                                                meson);

        return propagator;
    }

    gsl_complex calculateInverseNonDiagonalPropagator(double zeroMomentum, double threeMomentum, double gamma)
    {   
        gsl_complex propagator = calculateNonDiagonalPropagator(zeroMomentum, threeMomentum, gamma);

        return ( gsl_complex_inverse( propagator ) );
    }

    ComplexSquareMatrixGSL calculateDiagonalPropagator(double zeroMomentum, double threeMomentum, double gamma)
    {   
        ComplexSquareMatrixGSL propagator = 
        diagonalMesonPropagator(parametersNJL, temperature, 
                                upQuarkEffectiveChemicalPotential, downQuarkEffectiveChemicalPotential, strangeQuarkEffectiveChemicalPotential, 
                                upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                                zeroMomentum, threeMomentum, gamma, integralPrecision,
                                meson);

        return propagator;
    }

    vector<gsl_complex> calculateInverseDiagonalPropagatorEigenvalues(double zeroMomentum, double threeMomentum, double gamma)
    {   
        ComplexSquareMatrixGSL propagator = calculateDiagonalPropagator(zeroMomentum, threeMomentum, gamma);

        vector<gsl_complex> eigenvalues = calculateEigenvalues3By3ComplexMatrix( propagator.inverse() );

        return eigenvalues;
    }

    void calculateMesonMassAndWidth(double precision, MultiRootFindingMethod method, double mesonMassGuess, double mesonWidthGuess)
    {   
        double x[2];
        x[0] = mesonMassGuess + 1E-6; 
        x[1] = mesonWidthGuess + 1E-6;
        
        multiDimensionalRootFind(2, precision, &x[0], this, &SU3NJL3DCutoffMesonMassEquations, method);

        mesonMass = x[0];
        mesonWidth = x[1];
    }
};

double mesonStateMassAtMeltingPoint(double , double , double , mesonState );


#endif
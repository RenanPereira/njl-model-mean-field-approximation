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


int SU3NJL3DCutoffNonDiagonalMesonMassEquations(const gsl_vector *, void *, gsl_vector *);


class SU3NJL3DCutoffNonDiagonalMeson
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

    double mesonMass = 0.0/0.0;
    double mesonWidth = 0.0/0.0;

    vector<double> diagonalMesonMasses = {0.0/0.0, 0.0/0.0, 0.0/0.0};
    vector<double> diagonalMesonWidths = {0.0/0.0, 0.0/0.0, 0.0/0.0};

public:
    SU3NJL3DCutoffNonDiagonalMeson(SU3NJL3DCutoffParameters parametersNJLAux, double TAux, 
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

    SU3NJL3DCutoffNonDiagonalMeson(void* auxiliar)
    {   
        parametersNJL = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->parametersNJL;
        temperature = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->temperature;
        upQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->upQuarkEffectiveChemicalPotential;
        downQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->downQuarkEffectiveChemicalPotential;
        strangeQuarkEffectiveChemicalPotential = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->strangeQuarkEffectiveChemicalPotential;
        upQuarkEffectiveMass = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->upQuarkEffectiveMass;
        downQuarkEffectiveMass = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->downQuarkEffectiveMass;
        strangeQuarkEffectiveMass = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->strangeQuarkEffectiveMass;
        integralPrecision = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->integralPrecision;
        meson = ((class SU3NJL3DCutoffNonDiagonalMeson *)(auxiliar))->meson;
    };


    //non-diagonal mesons
    double getMesonMass(){ return mesonMass; }
    double getMesonWidth(){ return mesonWidth; }

    gsl_complex getPropagator(double zeroMomentum, double threeMomentum, double gamma)
    {   
        gsl_complex propagator;
        propagator = nonDiagonalMesonPropagator(parametersNJL, temperature, 
                                                upQuarkEffectiveChemicalPotential, downQuarkEffectiveChemicalPotential, strangeQuarkEffectiveChemicalPotential, 
                                                upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                                                zeroMomentum, threeMomentum, gamma, integralPrecision,
                                                meson);

        return propagator;
    }

    gsl_complex getInversePropagator(double zeroMomentum, double threeMomentum, double gamma)
    {   
        gsl_complex propagator = getPropagator(zeroMomentum, threeMomentum, gamma);

        return ( gsl_complex_inverse( propagator ) );
    }

    void calculateMesonMassAndWidth(double precision, MultiRootFindingMethod method, double mesonMassGuess, double mesonWidthGuess)
    {   
        double x[2];
        x[0] = mesonMassGuess; 
        x[1] = mesonWidthGuess; 
        
        multiDimensionalRootFind(2, precision, &x[0], this, &SU3NJL3DCutoffNonDiagonalMesonMassEquations, method);

        mesonMass = x[0];
        mesonWidth = x[1];
    }
};




#endif
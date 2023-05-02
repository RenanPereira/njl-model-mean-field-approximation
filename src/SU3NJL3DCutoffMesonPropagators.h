#ifndef SU3NJL3DCUTOFFMESONPROPAGATORS_H
#define SU3NJL3DCUTOFFMESONPROPAGATORS_H

#include <gsl/gsl_complex.h>
#include "ComplexSquareMatrixGSL.h"
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

#endif
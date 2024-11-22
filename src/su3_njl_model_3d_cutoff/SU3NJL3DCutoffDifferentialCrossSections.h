#ifndef SU3NJL3DCUTOFFDIFFERENTIALCROSSSECTIONS_H
#define SU3NJL3DCUTOFFDIFFERENTIALCROSSSECTIONS_H

#include <gsl/gsl_complex.h>
#include "SU3NJL3DCutoff.h"

using namespace std;


enum scatteringProcess { UDUD, DUDU, USUS, SUSU, DSDS, SDSD, UUUU, DDDD, SSSS,
                         UDBarUDBar, USBarUSBar, DSBarDSBar,
                         DUBarDUBar, SUBarSUBar, SDBarSDBar,
                         UUBarUUBar, UUBarDDBar, UUBarSSBar,
                         DDBarUUBar, DDBarDDBar, DDBarSSBar,
                         SSBarUUBar, SSBarDDBar, SSBarSSBar,
                         UBarUBarUBarUBar, DBarDBarDBarDBar, SBarSBarSBarSBar,
                         UBarDBarUBarDBar, UBarSBarUBarSBar };

string toString(scatteringProcess );

void inOutMassesGivenScatteringProcess(double , double , double , scatteringProcess , double& , double& , double& , double& );

void inChemicalPotentialsGivenScatteringProcess(double , double , double , scatteringProcess , double& , double& );

void outChemicalPotentialsGivenScatteringProcess(double , double , double , scatteringProcess , double& , double& );

double xij(double , double , double , double );

double sChannelZeroMomentum(double );

double sChannelThreeMomentum();

double tChannelZeroMomentum(double , double , double , double , double );

double tChannelThreeMomentum(double , double , double , double , double , double );

double uChannel(double , double , double , double , double , double );

double uChannelZeroMomentum(double , double , double , double , double );

double uChannelThreeMomentum(double , double , double , double , double , double );

double lambdaTriangle(double , double , double );

double cosScatteringAngle(double , double , double , double , double , double );

double sinScatteringAngle(double , double , double , double , double , double );

double dsigmadtNJLQuarkQuarkScattering(double , double , double , double , double , double , double , gsl_complex , gsl_complex , gsl_complex , gsl_complex );

double dsigmadtNJLQuarkAntiquarkScattering(double , double , double , double , double , double , double , gsl_complex , gsl_complex , gsl_complex , gsl_complex );

double dsigmadtUDUD(SU3NJL3DCutoffParameters , double , 
                    double , double , double , 
                    double , double , double , 
                    double , double , double );

double dsigmadtDUDU(SU3NJL3DCutoffParameters , double , 
                    double , double , double , 
                    double , double , double , 
                    double , double , double );

double dsigmadtUSUS(SU3NJL3DCutoffParameters , double , 
                    double , double , double , 
                    double , double , double , 
                    double , double , double );

double dsigmadtSUSU(SU3NJL3DCutoffParameters , double , 
                    double , double , double , 
                    double , double , double , 
                    double , double , double );

double dsigmadtDSDS(SU3NJL3DCutoffParameters , double , 
                    double , double , double , 
                    double , double , double , 
                    double , double , double );

double dsigmadtSDSD(SU3NJL3DCutoffParameters , double , 
                    double , double , double , 
                    double , double , double , 
                    double , double , double );

double dsigmadtUUUU(SU3NJL3DCutoffParameters , double , 
                    double , double , double , 
                    double , double , double , 
                    double , double , double );

double dsigmadtDDDD(SU3NJL3DCutoffParameters , double , 
                    double , double , double , 
                    double , double , double , 
                    double , double , double );

double dsigmadtSSSS(SU3NJL3DCutoffParameters , double , 
                    double , double , double , 
                    double , double , double , 
                    double , double , double );

double dsigmadtUDBarUDBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtDUBarDUBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtUSBarUSBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtSUBarSUBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtDSBarDSBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtSDBarSDBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtUUBarUUBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtUUBarDDBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtUUBarSSBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtDDBarUUBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtDDBarDDBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtDDBarSSBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtSSBarUUBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtSSBarDDBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double dsigmadtSSBarSSBar(SU3NJL3DCutoffParameters , double , 
                          double , double , double , 
                          double , double , double , 
                          double , double , double );

double differentialCrossSectionProcess12To34(SU3NJL3DCutoffParameters , double , 
                                             double , double , double , 
                                             double , double , double , 
                                             double , double , double , scatteringProcess );



#endif
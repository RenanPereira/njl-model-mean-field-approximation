#ifndef GENERALPHYSICSANDMATH_H
#define GENERALPHYSICSANDMATH_H


const double hc_MeVfm = 197.3269804;
const double hc_GeVfm = hc_MeVfm/1000.0;

const double electronMass_MeV = 0.511;
const double electronMass_GeV = electronMass_MeV/1000.0;


enum NJL3DCutoffRegularizationScheme { cutoffEverywhere, 
									   cutoffEverywhereWithCTmu,
									   cutoffOnDivergentIntegralsOnly };



double boseDistribution(double , double );

double fermiDistribution(double , double );

double fermiMomentum(double , double );

double heavisideTheta(double );

double puiseuxExpansionTln1plusExpArgOverT(double , double );

double sign(double );

#endif

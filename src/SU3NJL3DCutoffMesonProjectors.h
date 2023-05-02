#ifndef SU3NJL3DCUTOFFMESONPROJECTORS_H
#define SU3NJL3DCUTOFFMESONPROJECTORS_H

#include "SU3NJL3DCutoff.h"


class SU3NJL3DCutoffMesonProjector
{
private:
    double gS = 0.0;
	double kappa = 0.0;
	double g1 = 0.0;
	double g2 = 0.0;

    double sigmaU = 0.0/0.0;
    double sigmaD = 0.0/0.0;
    double sigmaS = 0.0/0.0;

public:
    SU3NJL3DCutoffMesonProjector(SU3NJL3DCutoffParameters parametersNJL, double T, double effChemPotU, double effChemPotD, double effChemPotS, double effMassU, double effMassD, double effMassS, double integralPrecision)
    {
        //necessary parameters
        NJLDimensionfulCouplings couplings = parametersNJL.getDimensionfulCouplings();

        NJL3DCutoffRegularizationScheme reguChoice = parametersNJL.getNJL3DCutoffRegularizationScheme();
        double cutoff = parametersNJL.getThreeMomentumCutoff();

        double Nc = parametersNJL.getNumberOfColours();

        //scalar-pseudoscalar couplings
        gS = couplings.getFourQuarkSPCoupling();
        kappa = couplings.getDeterminantCoupling();
        g1 = couplings.getEightQuarkSPOziViolatingCoupling();
        g2 = couplings.getEightQuarkSPNonOziViolatingCoupling();

        //Calculate the sigma field and density for each quark flavour
        sigmaU = sigmaNJL3DCutoff(reguChoice, cutoff, Nc, T, effChemPotU, effMassU, integralPrecision);
        sigmaD = sigmaNJL3DCutoff(reguChoice, cutoff, Nc, T, effChemPotD, effMassD, integralPrecision);
        sigmaS = sigmaNJL3DCutoff(reguChoice, cutoff, Nc, T, effChemPotS, effMassS, integralPrecision);
    };

    double pseudoscalar00();
    double pseudoscalar08();
    double pseudoscalar03();
    double pseudoscalar33();
    double pseudoscalar38();
    double pseudoscalar88();
    double pseudoscalar11();
    double pseudoscalar22();
    double pseudoscalar44();
    double pseudoscalar55();
    double pseudoscalar66();
    double pseudoscalar77();

    double scalar00();
    double scalar08();
    double scalar03();
    double scalar38();
    double scalar33();
    double scalar88();
    double scalar11();
    double scalar22();
    double scalar44();
    double scalar55();
    double scalar66();
    double scalar77();

};





#endif
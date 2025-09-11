#ifndef SU3NJL3DCUTOFFINTEGRATEDCROSSSECTIONS_H
#define SU3NJL3DCUTOFFINTEGRATEDCROSSSECTIONS_H

#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffDifferentialCrossSections.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCrossSections.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"

using namespace std;


class IntegratedCrossSectionIntegrand : public GeneralIntegrandParameters
{
private:
    string integralID = "notDefined";
    SU3NJL3DCutoffParameters parametersNJL;
    double temperature = 0.0/0.0;
    double upQuarkEffectiveChemicalPotential = 0.0/0.0;
    double downQuarkEffectiveChemicalPotential = 0.0/0.0;
    double strangeQuarkEffectiveChemicalPotential = 0.0/0.0;
    double upQuarkEffectiveMass = 0.0/0.0;
    double downQuarkEffectiveMass = 0.0/0.0;
    double strangeQuarkEffectiveMass = 0.0/0.0;
    double propagatorIntegralPrecision = 1E-8;
    scatteringProcess process;
    bool largeAngleScatteringContribution = false;
    double crossSectionIntegralPrecision = 0.0/0.0;
    double integratedCrossSectionIntegralPrecision_dXdYdZ = 0.0/0.0;
    double integratedCrossSectionIntegralPrecision_dXdY = 0.0/0.0;
    double centerOfMassEnergy = 0.0/0.0;
    double energy = 0.0/0.0;
    double momentumParticle1 = 0.0/0.0;
    double momentumParticle2 = 0.0/0.0;
    double normalizationRiemannSum_ds = 1.0;

public:
    IntegratedCrossSectionIntegrand(
        string integralIDAux, SU3NJL3DCutoffParameters parametersNJLAux, double temperatureAux, 
        double upQuarkEffectiveChemicalPotentialAux, double downQuarkEffectiveChemicalPotentialAux, double strangeQuarkEffectiveChemicalPotentialAux, 
        double upQuarkEffectiveMassAux, double downQuarkEffectiveMassAux, double strangeQuarkEffectiveMassAux, 
        double propagatorIntegralPrecisionAux, scatteringProcess processAux, 
        bool largeAngleScatteringContributionAux, double crossSectionIntegralPrecisionAux,
        double integratedCrossSectionIntegralPrecision_dXdYAux
    )
    {   
        integralID = integralIDAux;
        parametersNJL = parametersNJLAux;
        temperature = temperatureAux;
        upQuarkEffectiveChemicalPotential = upQuarkEffectiveChemicalPotentialAux;
        downQuarkEffectiveChemicalPotential = downQuarkEffectiveChemicalPotentialAux;
        strangeQuarkEffectiveChemicalPotential = strangeQuarkEffectiveChemicalPotentialAux;
        upQuarkEffectiveMass = upQuarkEffectiveMassAux;
        downQuarkEffectiveMass = downQuarkEffectiveMassAux;
        strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux;
        propagatorIntegralPrecision = propagatorIntegralPrecisionAux;
        process = processAux;
        largeAngleScatteringContribution = largeAngleScatteringContributionAux;
        crossSectionIntegralPrecision = crossSectionIntegralPrecisionAux;
        integratedCrossSectionIntegralPrecision_dXdY = integratedCrossSectionIntegralPrecision_dXdYAux;
    };

    IntegratedCrossSectionIntegrand(
    	string integralIDAux, SU3NJL3DCutoffParameters parametersNJLAux, double temperatureAux, 
        double upQuarkEffectiveChemicalPotentialAux, double downQuarkEffectiveChemicalPotentialAux, double strangeQuarkEffectiveChemicalPotentialAux, 
        double upQuarkEffectiveMassAux, double downQuarkEffectiveMassAux, double strangeQuarkEffectiveMassAux, 
        double propagatorIntegralPrecisionAux, scatteringProcess processAux, 
        bool largeAngleScatteringContributionAux, double crossSectionIntegralPrecisionAux,
        double integratedCrossSectionIntegralPrecision_dXdYdZAux, double integratedCrossSectionIntegralPrecision_dXdYAux
    )
    {   
        integralID = integralIDAux;
        parametersNJL = parametersNJLAux;
        temperature = temperatureAux;
        upQuarkEffectiveChemicalPotential = upQuarkEffectiveChemicalPotentialAux;
        downQuarkEffectiveChemicalPotential = downQuarkEffectiveChemicalPotentialAux;
        strangeQuarkEffectiveChemicalPotential = strangeQuarkEffectiveChemicalPotentialAux;
        upQuarkEffectiveMass = upQuarkEffectiveMassAux;
        downQuarkEffectiveMass = downQuarkEffectiveMassAux;
        strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux;
        propagatorIntegralPrecision = propagatorIntegralPrecisionAux;
        process = processAux;
        largeAngleScatteringContribution = largeAngleScatteringContributionAux;
        crossSectionIntegralPrecision = crossSectionIntegralPrecisionAux;
        integratedCrossSectionIntegralPrecision_dXdYdZ = integratedCrossSectionIntegralPrecision_dXdYdZAux;
        integratedCrossSectionIntegralPrecision_dXdY = integratedCrossSectionIntegralPrecision_dXdYAux;
    };

    IntegratedCrossSectionIntegrand(void* auxiliar)
    {   
        integralID = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->integralID;
        parametersNJL = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->parametersNJL;
        temperature = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->temperature;
        upQuarkEffectiveChemicalPotential = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->upQuarkEffectiveChemicalPotential;
        downQuarkEffectiveChemicalPotential = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->downQuarkEffectiveChemicalPotential;
        strangeQuarkEffectiveChemicalPotential = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->strangeQuarkEffectiveChemicalPotential;
        upQuarkEffectiveMass = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->upQuarkEffectiveMass;
        downQuarkEffectiveMass = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->downQuarkEffectiveMass;
        strangeQuarkEffectiveMass = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->strangeQuarkEffectiveMass;
        propagatorIntegralPrecision = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->propagatorIntegralPrecision;
        process = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->process;
        largeAngleScatteringContribution = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->largeAngleScatteringContribution;
        crossSectionIntegralPrecision = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->crossSectionIntegralPrecision;
        integratedCrossSectionIntegralPrecision_dXdYdZ = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->integratedCrossSectionIntegralPrecision_dXdYdZ;
        integratedCrossSectionIntegralPrecision_dXdY = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->integratedCrossSectionIntegralPrecision_dXdY;
        centerOfMassEnergy = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->centerOfMassEnergy;
        energy = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->energy;
        momentumParticle1 = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->momentumParticle1;
        momentumParticle2 = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->momentumParticle2;
        normalizationRiemannSum_ds = ((class IntegratedCrossSectionIntegrand *)(auxiliar))->normalizationRiemannSum_ds;
    };

    string getIntegralID(){ return integralID; }
    SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };
    double getTemperature(){ return temperature; };
    double getUpQuarkEffectiveMass(){ return upQuarkEffectiveMass; };
    double getDownQuarkEffectiveMass(){ return downQuarkEffectiveMass; };
    double getStrangeQuarkEffectiveMass(){ return strangeQuarkEffectiveMass; };
    double getUpQuarkEffectiveChemicalPotential(){ return upQuarkEffectiveChemicalPotential; };
    double getDownQuarkEffectiveChemicalPotential(){ return downQuarkEffectiveChemicalPotential; };
    double getStrangeQuarkEffectiveChemicalPotential(){ return strangeQuarkEffectiveChemicalPotential; };
    double getPropagatorIntegralPrecision(){ return propagatorIntegralPrecision; }
    scatteringProcess getProcess(){ return process; }
    bool getLargeAngleScatteringContribution(){ return largeAngleScatteringContribution; }
    double getCrossSectionIntegralPrecision(){ return crossSectionIntegralPrecision; }
    double getIntegratedCrossSectionIntegralPrecision_dXdYdZ(){ return integratedCrossSectionIntegralPrecision_dXdYdZ; }
    double getIntegratedCrossSectionIntegralPrecision_dXdY(){ return integratedCrossSectionIntegralPrecision_dXdY; }
    double getCenterOfMassEnergy(){ return centerOfMassEnergy; }
    double getMomentumParticle1(){ return momentumParticle1; }
    double getMomentumParticle2(){ return momentumParticle2; }
    double getNormalizationRiemannSum_ds(){ return normalizationRiemannSum_ds; }

    void setIntegralID(string integralIDAux){ integralID = integralIDAux; }
    void setCenterOfMassEnergy(double centerOfMassEnergyAux){ centerOfMassEnergy = centerOfMassEnergyAux; }
    void setEnergy(double energyAux){ energy = energyAux; }
    void setMomentumParticle1(double momentumParticle1Aux){ momentumParticle1 = momentumParticle1Aux; }
    void setMomentumParticle2(double momentumParticle2Aux){ momentumParticle2 = momentumParticle2Aux; }
    void setNormalizationRiemannSum_ds(double normalizationRiemannSum_dsAux){ normalizationRiemannSum_ds = normalizationRiemannSum_dsAux; }

    void printIntegrandVariables() override
    {   
        cout << "integralID = " << integralID << "\n";

        cout << "SU3NJL3DCutoffParameters are not being printed!" << "\n";

        cout << "T = " << temperature << "\n";
        cout << "effChemPotU = " << upQuarkEffectiveChemicalPotential << "\n";
        cout << "effChemPotD = " << downQuarkEffectiveChemicalPotential << "\n";
        cout << "effChemPotS = " << strangeQuarkEffectiveChemicalPotential << "\n";
        cout << "effMassU = " << upQuarkEffectiveMass << "\n";
        cout << "effMassD = " << downQuarkEffectiveMass << "\n";
        cout << "effMassS = " << strangeQuarkEffectiveMass << "\n";
        cout << "propagatorPrecision = " << propagatorIntegralPrecision << "\n";
        cout << "scatteringProcess = " << toString(process) << "\n";
        cout << "s = " << centerOfMassEnergy << "\n";
        cout << "E = " << energy << "\n";
        cout << "p1 = " << momentumParticle1 << "\n";
        cout << "p2 = " << momentumParticle2 << "\n";
    }
};

double fermiDiracIntegralIntegrand(double , void *);

double fermiDiracIntegral(double , double , double , double , double , double );

double fermiDiracIntegral(double , double , double , double , double );

double dsdEdepsilon_epsilonMinPlus(double , double , double , double );

double dsdEdepsilon_epsilonMinMinus(double , double , double , double );

double dsdEdepsilon_epsilonLambdaM1Plus(double , double , double );

double dsdEdepsilon_epsilonLambdaM1Minus(double , double , double );

double dsdEdepsilon_epsilonLambdaM2Plus(double , double , double );

double dsdEdepsilon_epsilonLambdaM2Minus(double , double , double );

double dsdEdepsilon_gamma1E(double , double , double , double );

double dsdEdepsilon_beta1E(double , double , double , double );

double dsdEdepsilon_gamma1epsilon(double , double , double , double );

double dsdEdepsilon_beta1epsilon(double , double , double , double );

double dsdEdepsilon_gamma2E(double , double , double , double );

double dsdEdepsilon_beta2E(double , double , double , double );

double dsdEdepsilon_gamma2epsilon(double , double , double , double );

double dsdEdepsilon_beta2epsilon(double , double , double , double );

double dsdEdepsilon_alpha1E(double );

double dsdEdepsilon_alpha1epsilon(double , double , double );

double dsdEdepsilon_ELambdaSwitchE(double , double , double );

double dsdEdepsilon_ELambdaswitchepsilon(double , double , double );

double dsdEdepsilon_epsilonMax(double , double , double , double , double );

double dsdEdepsilon_epsilonMin(double , double , double , double , double );

double dsdEdepsilon_EMinM1LargerM2(double , double , double , double );

double dsdEdepsilon_EMinM2LargerM1(double , double , double , double );

double dsdEdepsilon_EMin(double , double , double , double );

double dsdEdepsilon_EMaxM1LargerM2(double , double , double , double );

double dsdEdepsilon_EMaxM2LargerM1(double , double , double , double );

double dsdEdepsilon_EMax(double , double , double , double );

double integratedCrossSectionCOVVolumeIntegrand_dsdE(double , void *);

double integratedCrossSectionCOVVolumeIntegrand_ds(double , void *);

double integratedCrossSectionCOVVolume(double , double , double , double , double );

double integratedCrossSectionOGVolumeIntegrand_dp1dp2dtheta(double , void *);

double integratedCrossSectionOGVolumeIntegrand_dp1dp2(double , void *);

double integratedCrossSectionOGVolumeIntegrand_dp1(double , void *);

double integratedCrossSectionOGVolume(double , double , double , double , double , double );

double nEta(double , double , double , double , double , double , double , double );

double integratedCrossSectionIntegrand_dsdE(double , void *);

double integratedCrossSectionIntegrand_ds(double , void *);

double integratedCrossSectionCOVNormalizedIntegrand_dx(double , void *);

double integratedCrossSectionProcess12To34(
    SU3NJL3DCutoffParameters , double , 
    double , double , double , 
    double , double , double , 
    double , scatteringProcess , 
    bool , double , 
    double , double 
);

double integratedCrossSectionOGIntegrand_dp1dp2dtheta(double , void *);

double integratedCrossSectionOGIntegrand_dp1dp2(double , void *);

double integratedCrossSectionOGIntegrand_dp1(double , void *);

double integratedCrossSectionOGProcess12To34(
    SU3NJL3DCutoffParameters , double , 
    double , double , double , 
    double , double , double , 
    double , scatteringProcess , 
    bool , double , 
    double , double , double 
);

double probabilityKlevanskyHeavisideSolutionPlus(double , double , double , double );

double probabilityKlevanskyHeavisideSolutionMinus(double , double , double , double );

double probabilityKlevanskyIntegrand_dx(double , void *);

double probabilityKlevansky(
    SU3NJL3DCutoffParameters , double , 
    double , double , double , 
    double , double , double , 
    double , scatteringProcess ,
    double 
);

double probabilityKlevansky(double s, void *parameters);

double integratedCrossSectionKlevanskyIntegrand_ds(double , void *);

double integratedCrossSectionKlevanskyNormalizedIntegrand_dx(double , void *);

double integratedCrossSectionProcess12To34Klevansky(
    SU3NJL3DCutoffParameters , double , 
    double , double , double , 
    double , double , double , 
    double , scatteringProcess , 
    bool , double , 
    double , double 
);

double nonNormalizedProbabilityZhuangIntegrand_ds(double , void *);

double probabilityNormalizationInverseZhuang(
    SU3NJL3DCutoffParameters , double , 
    double , double , double , 
    double , double , double , 
    scatteringProcess ,
    double 
);

double integratedCrossSectionZhuangIntegrand_ds(double , void *);

double integratedCrossSectionZhuangNormalizedIntegrand_dx(double , void *);

double integratedCrossSectionProcess12To34Zhuang(
    SU3NJL3DCutoffParameters , double , 
    double , double , double , 
    double , double , double , 
    double , scatteringProcess , 
    bool , double , 
    double 
);


enum IntegratedCrossSectionApproximationMethod 
{ 
    COMPLETE_OG, 
    COMPLETE_COV, 
    KLEVANSKY, 
    ZHUANG 
};

inline const std::map<IntegratedCrossSectionApproximationMethod, std::string> IntegratedCrossSectionApproximationMethodMap = 
{
    {IntegratedCrossSectionApproximationMethod::COMPLETE_OG, "COMPLETE_OG"},
    {IntegratedCrossSectionApproximationMethod::COMPLETE_COV, "COMPLETE_COV"},
    {IntegratedCrossSectionApproximationMethod::KLEVANSKY, "KLEVANSKY"},
    {IntegratedCrossSectionApproximationMethod::ZHUANG, "ZHUANG"}
};


string toString(IntegratedCrossSectionApproximationMethod );

IntegratedCrossSectionApproximationMethod stringToIntegratedCrossSectionApproximationMethod(const std::string& );

bool isValidIntegratedCrossSectionApproximationMethod(const string& );

class SU3NJL3DCutoffIntegratedCrossSection
{
private:
    SU3NJL3DCutoffParameters parametersNJL;
    double temperature = 0.0/0.0;
    double upQuarkEffectiveChemicalPotential = 0.0/0.0;
    double downQuarkEffectiveChemicalPotential = 0.0/0.0;
    double strangeQuarkEffectiveChemicalPotential = 0.0/0.0;
    double upQuarkEffectiveMass = 0.0/0.0;
    double downQuarkEffectiveMass = 0.0/0.0;
    double strangeQuarkEffectiveMass = 0.0/0.0;
    double propagatorIntegralPrecision = 1E-8;
    scatteringProcess process;
    bool largeAngleScatteringContribution = false;
    double crossSectionIntegralPrecision = 1E-4;
    double integratedCrossSectionIntegralPrecision_dXdYdZ = 0.0/0.0;
    double integratedCrossSectionIntegralPrecision_dXdY = 1E-12;
    double integratedCrossSectionIntegralPrecision_dX = 1E-3;
    IntegratedCrossSectionApproximationMethod approximationMethod;

    double upQuarkNumber = 0.0/0.0;
    double downQuarkNumber = 0.0/0.0;
    double strangeQuarkNumber = 0.0/0.0;
    double upAntiquarkNumber = 0.0/0.0;
    double downAntiquarkNumber = 0.0/0.0;
    double strangeAntiquarkNumber = 0.0/0.0;

    double integratedCrossSection = 0.0/0.0; 
public:
    SU3NJL3DCutoffIntegratedCrossSection(){};
    SU3NJL3DCutoffIntegratedCrossSection(
        SU3NJL3DCutoffParameters parametersNJLAux, double temperatureAux, 
        double effChemPotUAux, double effChemPotDAux, double effChemPotSAux, 
        double effMassUAux, double effMassDAux, double effMassSAux, 
        double propagatorIntegralPrecisionAux, scatteringProcess processAux, 
        bool largeAngleScatteringContributionAux, double crossSectionIntegralPrecisionAux,
        double integratedCrossSectionIntegralPrecision_dXdYdZAux, 
        double integratedCrossSectionIntegralPrecision_dXdYAux, 
        double integratedCrossSectionIntegralPrecision_dXAux,
        IntegratedCrossSectionApproximationMethod approximationMethodAux
    )
    {
        parametersNJL = parametersNJLAux;
        temperature = temperatureAux;
        upQuarkEffectiveChemicalPotential = effChemPotUAux;
        downQuarkEffectiveChemicalPotential = effChemPotDAux;
        strangeQuarkEffectiveChemicalPotential = effChemPotSAux;
        upQuarkEffectiveMass = effMassUAux;
        downQuarkEffectiveMass = effMassDAux;
        strangeQuarkEffectiveMass = effMassSAux;
        propagatorIntegralPrecision = propagatorIntegralPrecisionAux;
        process = processAux;
        largeAngleScatteringContribution = largeAngleScatteringContributionAux;
        crossSectionIntegralPrecision = crossSectionIntegralPrecisionAux;
        integratedCrossSectionIntegralPrecision_dXdYdZ = integratedCrossSectionIntegralPrecision_dXdYdZAux;
        integratedCrossSectionIntegralPrecision_dXdY = integratedCrossSectionIntegralPrecision_dXdYAux;
        integratedCrossSectionIntegralPrecision_dX = integratedCrossSectionIntegralPrecision_dXAux;
        approximationMethod = approximationMethodAux;
        if ( approximationMethod!=COMPLETE_OG )
        {
            cout << "Calling a constructor for SU3NJL3DCutoffIntegratedCrossSection that is not appropriate for the chosen approximation method! Aborting!\n";
            abort();
        }
    }
    SU3NJL3DCutoffIntegratedCrossSection(
        SU3NJL3DCutoffParameters parametersNJLAux, double temperatureAux, 
        double effChemPotUAux, double effChemPotDAux, double effChemPotSAux, 
        double effMassUAux, double effMassDAux, double effMassSAux, 
        double propagatorIntegralPrecisionAux, scatteringProcess processAux, 
        bool largeAngleScatteringContributionAux, double crossSectionIntegralPrecisionAux,
        double integratedCrossSectionIntegralPrecision_dXdYAux, 
        double integratedCrossSectionIntegralPrecision_dXAux,
        IntegratedCrossSectionApproximationMethod approximationMethodAux
    )
    {
        parametersNJL = parametersNJLAux;
        temperature = temperatureAux;
        upQuarkEffectiveChemicalPotential = effChemPotUAux;
        downQuarkEffectiveChemicalPotential = effChemPotDAux;
        strangeQuarkEffectiveChemicalPotential = effChemPotSAux;
        upQuarkEffectiveMass = effMassUAux;
        downQuarkEffectiveMass = effMassDAux;
        strangeQuarkEffectiveMass = effMassSAux;
        propagatorIntegralPrecision = propagatorIntegralPrecisionAux;
        process = processAux;
        largeAngleScatteringContribution = largeAngleScatteringContributionAux;
        crossSectionIntegralPrecision = crossSectionIntegralPrecisionAux;
        integratedCrossSectionIntegralPrecision_dXdY = integratedCrossSectionIntegralPrecision_dXdYAux;
        integratedCrossSectionIntegralPrecision_dX = integratedCrossSectionIntegralPrecision_dXAux;
        approximationMethod = approximationMethodAux;
        if ( approximationMethod==COMPLETE_OG )
        {
            cout << "Calling a constructor for SU3NJL3DCutoffIntegratedCrossSection that is not appropriate for the chosen approximation method! Aborting!\n";
            abort();
        }
    }
    SU3NJL3DCutoffIntegratedCrossSection(
        SU3NJL3DCutoffParameters parametersNJLAux, double temperatureAux, 
        double effChemPotUAux, double effChemPotDAux, double effChemPotSAux, 
        double effMassUAux, double effMassDAux, double effMassSAux, 
        double propagatorIntegralPrecisionAux, scatteringProcess processAux, 
        bool largeAngleScatteringContributionAux, double crossSectionIntegralPrecisionAux,
        double integratedCrossSectionIntegralPrecision_dXAux,
        IntegratedCrossSectionApproximationMethod approximationMethodAux
    )
    {
        parametersNJL = parametersNJLAux;
        temperature = temperatureAux;
        upQuarkEffectiveChemicalPotential = effChemPotUAux;
        downQuarkEffectiveChemicalPotential = effChemPotDAux;
        strangeQuarkEffectiveChemicalPotential = effChemPotSAux;
        upQuarkEffectiveMass = effMassUAux;
        downQuarkEffectiveMass = effMassDAux;
        strangeQuarkEffectiveMass = effMassSAux;
        propagatorIntegralPrecision = propagatorIntegralPrecisionAux;
        process = processAux;
        largeAngleScatteringContribution = largeAngleScatteringContributionAux;
        crossSectionIntegralPrecision = crossSectionIntegralPrecisionAux;
        integratedCrossSectionIntegralPrecision_dX = integratedCrossSectionIntegralPrecision_dXAux;
        approximationMethod = approximationMethodAux;
        if ( approximationMethod!=ZHUANG )
        {
            cout << "Calling constructor for SU3NJL3DCutoffIntegratedCrossSection that is not appropriate for the chosen approximation method! Aborting!\n";
            abort();
        }
    }
    SU3NJL3DCutoffIntegratedCrossSection(
        SU3NJL3DCutoffParameters parametersNJLAux, double temperatureAux, 
        double effChemPotUAux, double effChemPotDAux, double effChemPotSAux, 
        double effMassUAux, double effMassDAux, double effMassSAux, 
        scatteringProcess processAux, 
        bool largeAngleScatteringContributionAux,
        IntegratedCrossSectionApproximationMethod approximationMethodAux
    )
    {
        parametersNJL = parametersNJLAux;
        temperature = temperatureAux;
        upQuarkEffectiveChemicalPotential = effChemPotUAux;
        downQuarkEffectiveChemicalPotential = effChemPotDAux;
        strangeQuarkEffectiveChemicalPotential = effChemPotSAux;
        upQuarkEffectiveMass = effMassUAux;
        downQuarkEffectiveMass = effMassDAux;
        strangeQuarkEffectiveMass = effMassSAux;
        process = processAux;
        largeAngleScatteringContribution = largeAngleScatteringContributionAux;
        approximationMethod = approximationMethodAux;

        if ( approximationMethod==COMPLETE_OG )
        {
            cout << "Calling a constructor for SU3NJL3DCutoffIntegratedCrossSection that is not appropriate for the chosen approximation method! Aborting!\n";
            abort();
        }
    }

    SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };
    double getTemperature(){ return temperature; };
    double getUpQuarkEffectiveMass(){ return upQuarkEffectiveMass; };
    double getDownQuarkEffectiveMass(){ return downQuarkEffectiveMass; };
    double getStrangeQuarkEffectiveMass(){ return strangeQuarkEffectiveMass; };
    double getUpQuarkEffectiveChemicalPotential(){ return upQuarkEffectiveChemicalPotential; };
    double getDownQuarkEffectiveChemicalPotential(){ return downQuarkEffectiveChemicalPotential; };
    double getStrangeQuarkEffectiveChemicalPotential(){ return strangeQuarkEffectiveChemicalPotential; };
    double getPropagatorIntegralPrecision(){ return propagatorIntegralPrecision; }
    scatteringProcess getProcess(){ return process; }
    bool getLargeAngleScatteringContribution(){ return largeAngleScatteringContribution; }
    double getCrossSectionIntegralPrecision(){ return crossSectionIntegralPrecision; }
    double getIntegratedCrossSectionIntegralPrecision_dXdYdZ(){ return integratedCrossSectionIntegralPrecision_dXdYdZ; }
    double getIntegratedCrossSectionIntegralPrecision_dXdY(){ return integratedCrossSectionIntegralPrecision_dXdY; }
    double getIntegratedCrossSectionIntegralPrecision_dX(){ return integratedCrossSectionIntegralPrecision_dX; }
    IntegratedCrossSectionApproximationMethod getApproximationMethod(){ return approximationMethod; }

    void setIntegratedCrossSection()
    {   
        if ( approximationMethod==COMPLETE_COV )
        {   
            integratedCrossSection = 
            integratedCrossSectionProcess12To34(
                parametersNJL, temperature, 
                upQuarkEffectiveChemicalPotential, 
                downQuarkEffectiveChemicalPotential, 
                strangeQuarkEffectiveChemicalPotential, 
                upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                propagatorIntegralPrecision, process, 
                largeAngleScatteringContribution, crossSectionIntegralPrecision,
                integratedCrossSectionIntegralPrecision_dXdY, integratedCrossSectionIntegralPrecision_dX
            );
        }
        else if( approximationMethod==KLEVANSKY )
        {
            integratedCrossSection = 
            integratedCrossSectionProcess12To34Klevansky(
                parametersNJL, temperature, 
                upQuarkEffectiveChemicalPotential, 
                downQuarkEffectiveChemicalPotential, 
                strangeQuarkEffectiveChemicalPotential, 
                upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                propagatorIntegralPrecision, process, 
                largeAngleScatteringContribution, crossSectionIntegralPrecision,
                integratedCrossSectionIntegralPrecision_dXdY, integratedCrossSectionIntegralPrecision_dX
            );
        }
        else if( approximationMethod==ZHUANG )
        {
            integratedCrossSection = 
            integratedCrossSectionProcess12To34Zhuang(
                parametersNJL, temperature, 
                upQuarkEffectiveChemicalPotential, 
                downQuarkEffectiveChemicalPotential, 
                strangeQuarkEffectiveChemicalPotential, 
                upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                propagatorIntegralPrecision, process, 
                largeAngleScatteringContribution, crossSectionIntegralPrecision,
                integratedCrossSectionIntegralPrecision_dX
            );
        }
        else if( approximationMethod==COMPLETE_OG )
        {
            integratedCrossSection = 
            integratedCrossSectionOGProcess12To34(
                parametersNJL, temperature, 
                upQuarkEffectiveChemicalPotential, 
                downQuarkEffectiveChemicalPotential, 
                strangeQuarkEffectiveChemicalPotential, 
                upQuarkEffectiveMass, downQuarkEffectiveMass, strangeQuarkEffectiveMass, 
                propagatorIntegralPrecision, process, 
                largeAngleScatteringContribution, crossSectionIntegralPrecision,
                integratedCrossSectionIntegralPrecision_dXdYdZ, 
                integratedCrossSectionIntegralPrecision_dXdY, 
                integratedCrossSectionIntegralPrecision_dX
            );
        }
    }

    void setQuarkNumbers()
    {   
        //get number of colors from parameters
        double Nc = parametersNJL.getNumberOfColours();
        
        //calculate quark numbers
        upQuarkNumber = fermiDiracIntegral(Nc, temperature, upQuarkEffectiveChemicalPotential, upQuarkEffectiveMass, propagatorIntegralPrecision);
        downQuarkNumber = fermiDiracIntegral(Nc, temperature, downQuarkEffectiveChemicalPotential, downQuarkEffectiveMass, propagatorIntegralPrecision);
        strangeQuarkNumber = fermiDiracIntegral(Nc, temperature, strangeQuarkEffectiveChemicalPotential, strangeQuarkEffectiveMass, propagatorIntegralPrecision);
        
        //calculate antiquark numbers
        upAntiquarkNumber = fermiDiracIntegral(Nc, temperature, -upQuarkEffectiveChemicalPotential, upQuarkEffectiveMass, propagatorIntegralPrecision);
        downAntiquarkNumber = fermiDiracIntegral(Nc, temperature, -downQuarkEffectiveChemicalPotential, downQuarkEffectiveMass, propagatorIntegralPrecision);
        strangeAntiquarkNumber = fermiDiracIntegral(Nc, temperature, -strangeQuarkEffectiveChemicalPotential, strangeQuarkEffectiveMass, propagatorIntegralPrecision);
    }

    double getUpQuarkNumber(){ return upQuarkNumber; }
    double getDownQuarkNumber(){ return downQuarkNumber; }
    double getStrangeQuarkNumber(){ return strangeQuarkNumber; }
    double getUpAntiquarkNumber(){ return upAntiquarkNumber; }
    double getDownAntiquarkNumber(){ return downAntiquarkNumber; }
    double getStrangeAntiquarkNumber(){ return strangeAntiquarkNumber; }
    double getIntegratedCrossSection(){ return integratedCrossSection; }
};


vector<SU3NJL3DCutoffIntegratedCrossSection> evaluateIntegratedCrossSectionAlongTrajectory(
    vector<SU3NJL3DCutoffFixedChemPotTemp> ,
    scatteringProcess , double , bool , 
    double , double , double , 
    IntegratedCrossSectionApproximationMethod 
);

void writeIntegratedCrossSectionToFile(vector<SU3NJL3DCutoffIntegratedCrossSection> , string );

void evaluateIntegratedCrossSectionAlongFixedChemicalPotentialTrajectory(
    vector<SU3NJL3DCutoffFixedChemPotTemp> ,
    scatteringProcess , double , bool , 
    double , double , double , 
    IntegratedCrossSectionApproximationMethod 
);

void evaluateAllIsospinSymmetricIntegratedCrossSectionsAlongFixedChemicalPotentialTrajectory(
    vector<SU3NJL3DCutoffFixedChemPotTemp> ,
    double , bool , 
    double , double , double , 
    IntegratedCrossSectionApproximationMethod 
);

void evaluateIntegratedCrossSectionAlongFixedTemperatureTrajectory(
    vector<SU3NJL3DCutoffFixedChemPotTemp> ,
    scatteringProcess , double , bool , 
    double , double , double , 
    IntegratedCrossSectionApproximationMethod 
);

void evaluateAllIsospinSymmetricIntegratedCrossSectionsAlongFixedTemperatureTrajectory(
    vector<SU3NJL3DCutoffFixedChemPotTemp> ,
    double , bool , 
    double , double , double , 
    IntegratedCrossSectionApproximationMethod 
);

void evaluateIsospinSymmetricIntegratedCrossSectionsWithZeroChemicalPotential(
    SU3NJL3DCutoffParameters& , 
	double , 
	MultiRootFindingMethod , 
	double , double , double , 
    MultiRootFindingMethod , 
    double , double , 
    int , int , 
    bool , 
    IntegratedCrossSectionApproximationMethod ,
    double , double , double , double 
);

void evaluateIntegratedCrossSectionsWithFixedTemperature(
    SU3NJL3DCutoffVacuum ,
    double , double , int , int , int , double , double , bool , 
    IntegratedCrossSectionApproximationMethod ,
    double , double , double , double 
);

void evaluateIntegratedCrossSectionsWithFixedChemicalPotential(
    SU3NJL3DCutoffVacuum ,
    double , double , double , double , int , int , int , bool , 
    IntegratedCrossSectionApproximationMethod ,
    double , double , double , double 
);

#endif
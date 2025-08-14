#ifndef SU3NJL3DCUTOFFCROSSSECTIONS_H
#define SU3NJL3DCUTOFFCROSSSECTIONS_H

#include "gsl_wrapper/root_solver_gsl.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoff.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffDifferentialCrossSections.h"

using namespace std;


class CrossSectionIntegrand : public GeneralIntegrandParameters
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
    double centerOfMassEnergy = 0.0/0.0;
    double propagatorIntegralPrecision = 1E-8;
    scatteringProcess process;
    bool largeAngleScatteringContribution = false;

public:
    CrossSectionIntegrand(string integralIDAux, SU3NJL3DCutoffParameters parametersNJLAux, double temperatureAux, 
                          double upQuarkEffectiveChemicalPotentialAux, double downQuarkEffectiveChemicalPotentialAux, double strangeQuarkEffectiveChemicalPotentialAux, 
                          double upQuarkEffectiveMassAux, double downQuarkEffectiveMassAux, double strangeQuarkEffectiveMassAux, 
                          double centerOfMassEnergyAux, double propagatorIntegralPrecisionAux, scatteringProcess processAux, 
                          bool largeAngleScatteringContributionAux)
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
        centerOfMassEnergy = centerOfMassEnergyAux;
        propagatorIntegralPrecision = propagatorIntegralPrecisionAux;
        process = processAux;
        largeAngleScatteringContribution = largeAngleScatteringContributionAux;
    };

    CrossSectionIntegrand(void* auxiliar)
    {   
        integralID = ((class CrossSectionIntegrand *)(auxiliar))->integralID;
        parametersNJL = ((class CrossSectionIntegrand *)(auxiliar))->parametersNJL;
        temperature = ((class CrossSectionIntegrand *)(auxiliar))->temperature;
        upQuarkEffectiveChemicalPotential = ((class CrossSectionIntegrand *)(auxiliar))->upQuarkEffectiveChemicalPotential;
        downQuarkEffectiveChemicalPotential = ((class CrossSectionIntegrand *)(auxiliar))->downQuarkEffectiveChemicalPotential;
        strangeQuarkEffectiveChemicalPotential = ((class CrossSectionIntegrand *)(auxiliar))->strangeQuarkEffectiveChemicalPotential;
        upQuarkEffectiveMass = ((class CrossSectionIntegrand *)(auxiliar))->upQuarkEffectiveMass;
        downQuarkEffectiveMass = ((class CrossSectionIntegrand *)(auxiliar))->downQuarkEffectiveMass;
        strangeQuarkEffectiveMass = ((class CrossSectionIntegrand *)(auxiliar))->strangeQuarkEffectiveMass;
        centerOfMassEnergy = ((class CrossSectionIntegrand *)(auxiliar))->centerOfMassEnergy;
        propagatorIntegralPrecision = ((class CrossSectionIntegrand *)(auxiliar))->propagatorIntegralPrecision;
        process = ((class CrossSectionIntegrand *)(auxiliar))->process;
        largeAngleScatteringContribution = ((class CrossSectionIntegrand *)(auxiliar))->largeAngleScatteringContribution;
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
    double getCenterOfMassEnergy(){ return centerOfMassEnergy; }
    double getPropagatorIntegralPrecision(){ return propagatorIntegralPrecision; }
    scatteringProcess getProcess(){ return process; }
    bool getLargeAngleScatteringContribution(){ return largeAngleScatteringContribution; }

    void setIntegralID(string integralIDAux){ integralID = integralIDAux; }

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
        cout << "s = " << centerOfMassEnergy << "\n";
        cout << "propagatorPrecision = " << propagatorIntegralPrecision << "\n";
        cout << "scatteringProcess = " << toString(process) << "\n";
    }
};


double momentumCM(double , double , double );

double energyOutParticle1(double , double , double );

double energyOutParticle2(double , double , double );

double energyOutParticle3(double , double , double );

double energyOutParticle4(double , double , double );

double relativeVelocityCM(double , double , double );

double relativeVelocity(double , double , double , double , double );

double tChannelMin(double , double , double , double , double );

double tChannelMax(double , double , double , double , double );

double centerOfMassEnergyThreshold(double , double , double , double );

double sMaximum(double , double , double , double , double );

double sMaximumKlevansky(double , double , double , double );

double crossSectionProcess12To34Integrand(double x, void *parameters);

double crossSectionProcess12To34(SU3NJL3DCutoffParameters , double , 
                                 double , double , double , 
                                 double , double , double , 
                                 double , double , scatteringProcess ,
                                 bool , double );

void evaluateCrossSectionProcess12To34ToFile(SU3NJL3DCutoffParameters , double , 
                                             double , double , double , 
                                             double , double , double , 
                                             double , scatteringProcess , 
                                             bool , double ,
                                             int , int );

void evaluateCrossSectionsKlevanskyPaper(SU3NJL3DCutoffParameters , double , 
                                         double , double , double , 
                                         double , double , double , 
                                         double , 
                                         bool , double ,
                                         int , int );

void evaluateCrossSectionsEqualLightMassesEqualChemicalPotential(
    SU3NJL3DCutoffParameters , double , 
    double , double , double , 
    double , double , double , 
    double , 
    bool , double ,
    int ,
    int
);



#endif
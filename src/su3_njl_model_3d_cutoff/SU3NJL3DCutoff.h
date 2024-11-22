#ifndef SU3NJL3DCUTOFF_H
#define SU3NJL3DCUTOFF_H

#include "Integration1DimGSL.h"
#include "NJLDimensionfulCouplings.h"

using namespace std;



class SU3NJL3DCutoffParameters
{
private:
    NJLDimensionfulCouplings couplings;

    NJL3DCutoffRegularizationScheme reguScheme;
    double threeMomentumCutoff = 0.0;

    double numberOfColours = 3;

    double sigmaIntegralPrecision = 1E-8;
    double thermoIntegralPrecision = 1E-8;

    double upQuarkCurrentMass = 0.0;
    double downQuarkCurrentMass = 0.0;
    double strangeQuarkCurrentMass = 0.0;

    string parameterSetName = "";

public:
    SU3NJL3DCutoffParameters(){};
    SU3NJL3DCutoffParameters(NJL3DCutoffRegularizationScheme reguSchemeAux, double threeMomentumCutoffAux, NJLDimensionfulCouplings couplingsAux, double m0u, double m0d, double m0s)
    {
        reguScheme = reguSchemeAux;
        threeMomentumCutoff = threeMomentumCutoffAux;
        couplings = couplingsAux;
        upQuarkCurrentMass = m0u;
        downQuarkCurrentMass = m0d;
        strangeQuarkCurrentMass = m0s;
    }
    SU3NJL3DCutoffParameters(double threeMomentumCutoffAux)
    {
        threeMomentumCutoff = threeMomentumCutoffAux;
    }

    void setSigmaIntegralPrecision(double sigmaIntegralPrecisionAux){ sigmaIntegralPrecision = sigmaIntegralPrecisionAux; };
    void setThermoIntegralPrecision(double thermoIntegralPrecisionAux){ thermoIntegralPrecision = thermoIntegralPrecisionAux; };
    void setNumberOfColours(double numberOfColoursAux){ numberOfColours = numberOfColoursAux; };
    void setParameterSetName(string parameterSetNameAux){ parameterSetName = parameterSetNameAux; };

    NJLDimensionfulCouplings getDimensionfulCouplings(){ return couplings; };
    NJL3DCutoffRegularizationScheme getNJL3DCutoffRegularizationScheme(){ return reguScheme; };
    double getThreeMomentumCutoff(){ return threeMomentumCutoff; };
    double getNumberOfColours(){ return numberOfColours; };
    double getSigmaIntegralPrecision(){ return sigmaIntegralPrecision; };
    double getThermoIntegralPrecision(){ return thermoIntegralPrecision; };
    double getUpQuarkCurrentMass(){ return upQuarkCurrentMass; };
    double getDownQuarkCurrentMass(){ return downQuarkCurrentMass; };
    double getStrangeQuarkCurrentMass(){ return strangeQuarkCurrentMass; };
    string getParameterSetName(){ return parameterSetName; };
};


class ThermodynamicsIntegrandParameters : public GeneralIntegrandParameters
{
public:
    string integralID;
    double temperature;
    double effectiveChemicalPotential;
    double chemicalPotential;
    double effectiveMass;

public:
    ThermodynamicsIntegrandParameters(string integralIDAux, double temperatureAux, double effectiveChemicalPotentialAux, double effectiveMassAux)
    {   
        integralID = integralIDAux;
        temperature = temperatureAux;
        effectiveChemicalPotential = effectiveChemicalPotentialAux;
        chemicalPotential = 0.0/0.0;
        effectiveMass = effectiveMassAux;
    };

    ThermodynamicsIntegrandParameters(string integralIDAux, double temperatureAux, double effectiveChemicalPotentialAux, double chemicalPotentialAux, double effectiveMassAux)
    {   
        integralID = integralIDAux;
        temperature = temperatureAux;
        effectiveChemicalPotential = effectiveChemicalPotentialAux;
        chemicalPotential = chemicalPotentialAux;
        effectiveMass = effectiveMassAux;
    };

    ThermodynamicsIntegrandParameters(void* auxiliar)
    {   
        integralID = ((class ThermodynamicsIntegrandParameters *)(auxiliar))->integralID;
        temperature = ((class ThermodynamicsIntegrandParameters *)(auxiliar))->temperature;
        effectiveChemicalPotential = ((class ThermodynamicsIntegrandParameters *)(auxiliar))->effectiveChemicalPotential;
        chemicalPotential = ((class ThermodynamicsIntegrandParameters *)(auxiliar))->chemicalPotential;
        effectiveMass = ((class ThermodynamicsIntegrandParameters *)(auxiliar))->effectiveMass;
    };

    double getTemperature(){ return temperature; };
    double getEffectiveChemicalPotential(){ return effectiveChemicalPotential; }
    double getChemicalPotential(){ return chemicalPotential; }
    double getEffectiveMass(){ return effectiveMass; };

    void printIntegrandVariables() override
    {   
        cout << integralID << "\n";
        cout << "T = " << temperature << "\n";
        cout << "effChemPot = " << effectiveChemicalPotential << "\n";
        cout << "chemPot = " << chemicalPotential << "\n";
        cout << "effMass = " << effectiveMass << "\n";
        cout << "(T, effChemPot, chemPot, effMass)\n";
        cout << "(" << temperature << ", " << effectiveChemicalPotential << ", " << chemicalPotential << ", " << effectiveMass << ")" << "\n";
    }
};


double SU3BaryonDensity(double , double , double );


double SU3NJLNulledGapEquation(NJLDimensionfulCouplings , double , double , double , double , double , double , double );

double SU3NJLQuarkChemicalPotential(NJLDimensionfulCouplings , double , double , double , double , double , double , double );

double SU3NJLInteractionPotential(NJLDimensionfulCouplings , double , double , double , double , double , double );


double fermionParticleDensityIntegrand(double , void *);

double fermionParticleDensity3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double );


double fermionPressureDivergentPrimitive3DCutoff(double , double );

double fermionPressureConvergentIntegrand(double , void *);

double fermionPressure3DCutoffStefanBoltzmannCTmu(double , double , double , double );

double fermionPressure3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double );


double fermionEnergyDensityConvergentIntegrand(double , void *);

double fermionEnergyDensity3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double , double );

double fermionEnergyDensity3DCutoffStefanBoltzmannCTmu(double , double , double , double , double );


double fermionEntropyDensityConvergentIntegrand(double , void *);

double fermionEntropyDensity3DCutoff(NJL3DCutoffRegularizationScheme , double , double , double , double , double );

double fermionEntropyDensity3DCutoffStefanBoltzmannCTmu(double , double , double , double );


double SU3NJL3DCutoffPressure(SU3NJL3DCutoffParameters , double , double , double , double , double , double , double );

double SU3NJL3DCutoffEnergyDensity(SU3NJL3DCutoffParameters , double , double , double , double , double , double , double );

double SU3NJL3DCutoffEntropyDensity(SU3NJL3DCutoffParameters , double , double , double , double , double , double , double );




#endif
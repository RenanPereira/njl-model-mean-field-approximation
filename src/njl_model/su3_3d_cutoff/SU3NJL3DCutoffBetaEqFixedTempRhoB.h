#ifndef SU3NJL3DCUTOFFBETAEQFIXEDTEMPRHOB_H
#define SU3NJL3DCUTOFFBETAEQFIXEDTEMPRHOB_H

#include <vector>
#include "generalPhysicsAndMath.h"
#include "njl_model/NJLDimensionfulCouplings.h"
using namespace std;



class SU3NJL3DCutoffBetaEqFixedTempRhoB
{
private:
	SU3NJL3DCutoffParameters parametersNJL;
	double electronMass;
	double temperature;
	double baryonDensity;

	double upQuarkEffectiveMass = 0.0/0.0;
	double downQuarkEffectiveMass = 0.0/0.0;
	double strangeQuarkEffectiveMass = 0.0/0.0;
	double upQuarkEffectiveChemicalPotential = 0.0/0.0;
	double downQuarkEffectiveChemicalPotential = 0.0/0.0;
	double strangeQuarkEffectiveChemicalPotential = 0.0/0.0;

	//Outputs
	double upQuarkSigma = 0.0/0.0;
	double downQuarkSigma = 0.0/0.0;
	double strangeQuarkSigma = 0.0/0.0;
	
	double upQuarkDensity = 0.0/0.0;
	double downQuarkDensity = 0.0;
	double strangeQuarkDensity = 0.0/0.0;

	double upQuarkChemicalPotential = 0.0/0.0;
	double downQuarkChemicalPotential = 0.0/0.0;
	double strangeQuarkChemicalPotential = 0.0/0.0;
	double baryonChemicalPotential = 0.0/0.0;

	double electronChemicalPotential = 0.0/0.0;
	double electronDensity = 0.0/0.0;

	//Thermodynamics
	double betaEqPressure = 0.0/0.0;
	double betaEqEnergyDensity = 0.0/0.0;
	double betaEqEntropyDensity = 0.0/0.0;

public:
	SU3NJL3DCutoffBetaEqFixedTempRhoB(void* );
	SU3NJL3DCutoffBetaEqFixedTempRhoB(SU3NJL3DCutoffParameters , double , double , double );
	SU3NJL3DCutoffBetaEqFixedTempRhoB(SU3NJL3DCutoffParameters , double , double , double , double , double , double , double , double );

	SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };

	double getElectronMass(){ return electronMass; };
	void setElectronMass(double electronMassAux){ electronMass = electronMassAux; };

	double getTemperature(){ return temperature; };
	double getBaryonDensity(){ return baryonDensity; };

	void setTemperature(double temperatureAux){ temperature = temperatureAux; };
	void setBaryonDensity(double baryonDensityAux){ baryonDensity = baryonDensityAux; };

	double getUpQuarkEffectiveMass(){ return upQuarkEffectiveMass; };
	double getDownQuarkEffectiveMass(){ return downQuarkEffectiveMass; };
	double getStrangeQuarkEffectiveMass(){ return strangeQuarkEffectiveMass; };
	double getUpQuarkEffectiveChemicalPotential(){ return upQuarkEffectiveChemicalPotential; };
	double getDownQuarkEffectiveChemicalPotential(){ return downQuarkEffectiveChemicalPotential; };
	double getStrangeQuarkEffectiveChemicalPotential(){ return strangeQuarkEffectiveChemicalPotential; };


	//solving the model
	void setSigmasAndDensitiesAndChemicalPotentials();
	void solve(double , MultiRootFindingMethod , double , double , double , double , double , double );
	bool testSolution(double );

	//get other quantities
	double getUpQuarkSigma(){ return upQuarkSigma; };
	double getDownQuarkSigma(){ return downQuarkSigma; };
	double getStrangeQuarkSigma(){ return strangeQuarkSigma; };

	double getUpQuarkDensity(){ return upQuarkDensity; };
	double getDownQuarkDensity(){ return downQuarkDensity; };
	double getStrangeQuarkDensity(){ return strangeQuarkDensity; };

	double getUpQuarkChemicalPotential(){ return upQuarkChemicalPotential; };
	double getDownQuarkChemicalPotential(){ return downQuarkChemicalPotential; };
	double getStrangeQuarkChemicalPotential(){ return strangeQuarkChemicalPotential; };

	double getElectronChemicalPotential(){ return electronChemicalPotential; };
	double getElectronDensity(){ return electronDensity; };

	//baryon chemical potential
	void setBaryonChemicalPotential(){ baryonChemicalPotential = getUpQuarkChemicalPotential() + 2*getDownQuarkChemicalPotential(); };
	double getBaryonChemicalPotential(){ return baryonChemicalPotential; }


	//thermodynamics
	double calculatePressureQuarks(double );
	double calculatePressureElectrons(double );
	double calculatePressure(double , double );
	void setBetaEqPressure(double , double );
	double getBetaEqPressure(){ return betaEqPressure; };

	double calculateEnergyDensityQuarks(double );
	double calculateEnergyDensityElectrons(double );
	double calculateEnergyDensity(double , double );
	void setBetaEqEnergyDensity(double , double );
	double getBetaEqEnergyDensity(){ return betaEqEnergyDensity; };

	double calculateEntropyDensityQuarks();
	double calculateEntropyDensityElectrons();
	double calculateEntropyDensity();
	void setBetaEqEntropyDensity();
	double getBetaEqEntropyDensity(){ return betaEqEntropyDensity; };

	void setBetaEqThermodynamics(double , double );


private:
	void setUpQuarkEffectiveMass(double upQuarkEffectiveMassAux){ upQuarkEffectiveMass = upQuarkEffectiveMassAux; };
	void setDownQuarkEffectiveMass(double downQuarkEffectiveMassAux){ downQuarkEffectiveMass = downQuarkEffectiveMassAux; };
	void setStrangeQuarkEffectiveMass(double strangeQuarkEffectiveMassAux){ strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux; };
	void setUpQuarkEffectiveChemicalPotential(double upQuarkEffectiveChemicalPotentialAux){ upQuarkEffectiveChemicalPotential = upQuarkEffectiveChemicalPotentialAux; };
	void setDownQuarkEffectiveChemicalPotential(double downQuarkEffectiveChemicalPotentialAux){ downQuarkEffectiveChemicalPotential = downQuarkEffectiveChemicalPotentialAux; };
	void setStrangeQuarkEffectiveChemicalPotential(double strangeQuarkEffectiveChemicalPotentialAux){ strangeQuarkEffectiveChemicalPotential = strangeQuarkEffectiveChemicalPotentialAux; };

	void setUpQuarkSigma(double upQuarkSigmaAux){ upQuarkSigma = upQuarkSigmaAux; };
	void setDownQuarkSigma(double downQuarkSigmaAux){ downQuarkSigma = downQuarkSigmaAux; };
	void setStrangeQuarkSigma(double strangeQuarkSigmaAux){ strangeQuarkSigma = strangeQuarkSigmaAux; };

	void setUpQuarkDensity(double upQuarkDensityAux){ upQuarkDensity = upQuarkDensityAux; };
	void setDownQuarkDensity(double downQuarkDensityAux){ downQuarkDensity = downQuarkDensityAux; };
	void setStrangeQuarkDensity(double strangeQuarkDensityAux){ strangeQuarkDensity = strangeQuarkDensityAux; };

	void setUpQuarkChemicalPotential(double upQuarkChemicalPotentialAux){ upQuarkChemicalPotential = upQuarkChemicalPotentialAux; };
	void setDownQuarkChemicalPotential(double downQuarkChemicalPotentialAux){ downQuarkChemicalPotential = downQuarkChemicalPotentialAux; };
	void setStrangeQuarkChemicalPotential(double strangeQuarkChemicalPotentialAux){ strangeQuarkChemicalPotential = strangeQuarkChemicalPotentialAux; };

	void setElectronChemicalPotential(double electronChemicalPotentialAux){ electronChemicalPotential = electronChemicalPotentialAux; };
	void setElectronDensity(double electronDensityAux){ electronDensity = electronDensityAux; };
};


int SU3NJL3DCutoffGapEquationsBetaEquilibriumFixedTemperature(const gsl_vector *, void *, gsl_vector *);

void writeSolutionsToFile(vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> , string , bool );

void writeEOSToFile(vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> , string , bool );

void writeEOSToFile(vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> , string , bool , double );

int SU3NJL3DCutoffChiralTransitionPointBetaEquilibriumFixedTemperature(const gsl_vector *, void *, gsl_vector *);

vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> findChiralTransitionPointsFixedTemperature(vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> , double , MultiRootFindingMethod );

void addVacuumSolution(SU3NJL3DCutoffVacuum , double , double , double , vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> &);

vector<SU3NJL3DCutoffBetaEqFixedTempRhoB> calculateZeroTemperatureSolutions(SU3NJL3DCutoffVacuum , double , double , int , double , MultiRootFindingMethod );

void writeBetaEquilibriumEOSAtZeroTemperatureToFile(SU3NJL3DCutoffVacuum , double , double , int , double , MultiRootFindingMethod , string );

#endif
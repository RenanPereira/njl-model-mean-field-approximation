#ifndef SU3NJL3DCutoffEqualChemPotFixedTempRhoB_H
#define SU3NJL3DCUTOFFEQUALCHEMPOTFIXEDTEMPRHOB_H

#include <vector>
#include "generalPhysicsAndMath.h"
#include "NJLDimensionfulCouplings.h"
using namespace std;


class SU3NJL3DCutoffEqualChemPotFixedTempRhoB
{
private:
	SU3NJL3DCutoffParameters parametersNJL;

	double temperature = 0.0;
	double baryonDensity = 0.0;

	double upQuarkEffectiveMass = 0.0/0.0;
	double downQuarkEffectiveMass = 0.0/0.0;
	double strangeQuarkEffectiveMass = 0.0/0.0;
	double quarkEffectiveChemicalPotential = 0.0/0.0;


public:
	SU3NJL3DCutoffEqualChemPotFixedTempRhoB(){};
	SU3NJL3DCutoffEqualChemPotFixedTempRhoB(void* );
	SU3NJL3DCutoffEqualChemPotFixedTempRhoB(SU3NJL3DCutoffParameters );

	SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };

	double getTemperature(){ return temperature; };
	double getBaryonDensity(){ return baryonDensity; };

	void setTemperature(double temperatureAux){ temperature = temperatureAux; };
	void setBaryonDensity(double baryonDensityAux){ baryonDensity = baryonDensityAux; };

	double getUpQuarkEffectiveMass(){ return upQuarkEffectiveMass; };
	double getDownQuarkEffectiveMass(){ return downQuarkEffectiveMass; };
	double getStrangeQuarkEffectiveMass(){ return strangeQuarkEffectiveMass; };

	double getQuarkEffectiveChemicalPotential(){ return quarkEffectiveChemicalPotential; };

	void solve(double , MultiRootFindingMethod , double , double , double , double );
	bool testSolution(double );

	double calculatePressure(double );
	double calculateEnergyDensity(double );
	double calculateEntropyDensity();

private:
	void setUpQuarkEffectiveMass(double upQuarkEffectiveMassAux){ upQuarkEffectiveMass = upQuarkEffectiveMassAux; };
	void setDownQuarkEffectiveMass(double downQuarkEffectiveMassAux){ downQuarkEffectiveMass = downQuarkEffectiveMassAux; };
	void setStrangeQuarkEffectiveMass(double strangeQuarkEffectiveMassAux){ strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux; };
	void setQuarkEffectiveChemicalPotential(double quarkEffectiveChemicalPotentialAux){ quarkEffectiveChemicalPotential = quarkEffectiveChemicalPotentialAux; };
};


int SU3NJL3DCutoffGapEquationsEqualChemicalPotentialFixedTemperatureBaryonDensity(const gsl_vector *, void *, gsl_vector *);



#endif
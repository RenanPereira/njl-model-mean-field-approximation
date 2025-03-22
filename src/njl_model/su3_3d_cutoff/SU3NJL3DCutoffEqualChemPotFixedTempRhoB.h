#ifndef SU3NJL3DCUTOFFEQUALCHEMPOTFIXEDTEMPRHOB_H
#define SU3NJL3DCUTOFFEQUALCHEMPOTFIXEDTEMPRHOB_H

#include <vector>
#include "gsl_wrapper/root_solver_gsl.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"


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

	double pressure = 0.0/0.0;
	double energyDensity = 0.0/0.0;
	double entropyDensity = 0.0/0.0;

public:
	SU3NJL3DCutoffEqualChemPotFixedTempRhoB(){};
	SU3NJL3DCutoffEqualChemPotFixedTempRhoB(void* );
	SU3NJL3DCutoffEqualChemPotFixedTempRhoB(SU3NJL3DCutoffParameters );
	SU3NJL3DCutoffEqualChemPotFixedTempRhoB(SU3NJL3DCutoffVacuum );

	SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };

	double getTemperature(){ return temperature; };
	double getBaryonDensity(){ return baryonDensity; };

	void setTemperature(double temperatureAux){ temperature = temperatureAux; };
	void setBaryonDensity(double baryonDensityAux){ baryonDensity = baryonDensityAux; };
	void setPressure(double pressureAux){ pressure = pressureAux; };
	void setEnergyDensity(double energyDensityAux){ energyDensity = energyDensityAux; };
	void setEntropyDensity(double entropyDensityAux){ entropyDensity = entropyDensityAux; };

	double getUpQuarkEffectiveMass(){ return upQuarkEffectiveMass; };
	double getDownQuarkEffectiveMass(){ return downQuarkEffectiveMass; };
	double getStrangeQuarkEffectiveMass(){ return strangeQuarkEffectiveMass; };
	double getPressure(){ return pressure; };
	double getEnergyDensity(){ return energyDensity; };
	double getEntropyDensity(){ return entropyDensity; };

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


vector<SU3NJL3DCutoffEqualChemPotFixedTempRhoB> 
solveFromVacuumToFiniteBaryonDensity(SU3NJL3DCutoffVacuum , 
                                     double , double , int , 
                                     double , MultiRootFindingMethod );

void writeSolutionsToFile(vector<SU3NJL3DCutoffEqualChemPotFixedTempRhoB> , string , bool );

#endif
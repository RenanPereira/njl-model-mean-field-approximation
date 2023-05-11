#ifndef SU3NJL3DCUTOFFFIXEDCHEMPOTTEMP_H
#define SU3NJL3DCUTOFFFIXEDCHEMPOTTEMP_H

#include <vector>
#include "generalPhysicsAndMath.h"
#include "NJLDimensionfulCouplings.h"
#include "SU3NJL3DCutoffVacuum.h"

using namespace std;


class SU3NJL3DCutoffFixedChemPotTemp
{
private:
	SU3NJL3DCutoffParameters parametersNJL;

	double temperature;
	double upQuarkChemicalPotential;
	double downQuarkChemicalPotential;
	double strangeQuarkChemicalPotential;

	double upQuarkEffectiveMass = 0.0/0.0;
	double downQuarkEffectiveMass = 0.0/0.0;
	double strangeQuarkEffectiveMass = 0.0/0.0;

public:
	SU3NJL3DCutoffFixedChemPotTemp(){};
	SU3NJL3DCutoffFixedChemPotTemp(void* );
	SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffParameters , double , double , double , double );

	SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };

	double getTemperature(){ return temperature; };
	void setTemperature(double temperatureAux){ temperature = temperatureAux; };

	double getUpQuarkEffectiveMass(){ return upQuarkEffectiveMass; };
	double getDownQuarkEffectiveMass(){ return downQuarkEffectiveMass; };
	double getStrangeQuarkEffectiveMass(){ return strangeQuarkEffectiveMass; };
	double getUpQuarkChemicalPotential(){ return upQuarkChemicalPotential; };
	double getDownQuarkChemicalPotential(){ return downQuarkChemicalPotential; };
	double getStrangeQuarkChemicalPotential(){ return strangeQuarkChemicalPotential; };

	void solve(double , MultiRootFindingMethod , double , double , double );
	bool testSolution(double );

private:
	void setUpQuarkEffectiveMass(double upQuarkEffectiveMassAux){ upQuarkEffectiveMass = upQuarkEffectiveMassAux; };
	void setDownQuarkEffectiveMass(double downQuarkEffectiveMassAux){ downQuarkEffectiveMass = downQuarkEffectiveMassAux; };
	void setStrangeQuarkEffectiveMass(double strangeQuarkEffectiveMassAux){ strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux; };
	void setUpQuarkChemicalPotential(double upQuarkChemicalPotentialAux){ upQuarkChemicalPotential = upQuarkChemicalPotentialAux; };
	void setDownQuarkChemicalPotential(double downQuarkChemicalPotentialAux){ downQuarkChemicalPotential = downQuarkChemicalPotentialAux; };
	void setStrangeQuarkChemicalPotential(double strangeQuarkChemicalPotentialAux){ strangeQuarkChemicalPotential = strangeQuarkChemicalPotentialAux; };
};


int SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature(const gsl_vector *, void *, gsl_vector *);

std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(SU3NJL3DCutoffVacuum , double , int , double , MultiRootFindingMethod );

std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromFiniteTemperatureToFiniteChemicalPotential(SU3NJL3DCutoffFixedChemPotTemp , double , int , double , MultiRootFindingMethod );

#endif
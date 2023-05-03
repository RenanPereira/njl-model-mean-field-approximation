#ifndef SU3NJL3DCUTOFFFIXEDCHEMPOTTEMP_H
#define SU3NJL3DCUTOFFFIXEDCHEMPOTTEMP_H

#include <vector>
#include "generalPhysicsAndMath.h"
#include "NJLDimensionfulCouplings.h"
using namespace std;


class SU3NJL3DCutoffFixedChemPotTemp
{
private:
	SU3NJL3DCutoffParameters parametersNJL;

	double temperature;
	double upQuarkEffectiveChemicalPotential;
	double downQuarkEffectiveChemicalPotential;
	double strangeQuarkEffectiveChemicalPotential;

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
	double getUpQuarkEffectiveChemicalPotential(){ return upQuarkEffectiveChemicalPotential; };
	double getDownQuarkEffectiveChemicalPotential(){ return downQuarkEffectiveChemicalPotential; };
	double getStrangeQuarkEffectiveChemicalPotential(){ return strangeQuarkEffectiveChemicalPotential; };

	void solve(double , MultiRootFindingMethod , double , double , double );
	bool testSolution(double );

private:
	void setUpQuarkEffectiveMass(double upQuarkEffectiveMassAux){ upQuarkEffectiveMass = upQuarkEffectiveMassAux; };
	void setDownQuarkEffectiveMass(double downQuarkEffectiveMassAux){ downQuarkEffectiveMass = downQuarkEffectiveMassAux; };
	void setStrangeQuarkEffectiveMass(double strangeQuarkEffectiveMassAux){ strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux; };
	void setUpQuarkEffectiveChemicalPotential(double upQuarkEffectiveChemicalPotentialAux){ upQuarkEffectiveChemicalPotential = upQuarkEffectiveChemicalPotentialAux; };
	void setDownQuarkEffectiveChemicalPotential(double downQuarkEffectiveChemicalPotentialAux){ downQuarkEffectiveChemicalPotential = downQuarkEffectiveChemicalPotentialAux; };
	void setStrangeQuarkEffectiveChemicalPotential(double strangeQuarkEffectiveChemicalPotentialAux){ strangeQuarkEffectiveChemicalPotential = strangeQuarkEffectiveChemicalPotentialAux; };
};


int SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature(const gsl_vector *, void *, gsl_vector *);


#endif
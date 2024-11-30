#ifndef SU3NJL3DCUTOFFVACUUM_H
#define SU3NJL3DCUTOFFVACUUM_H

#include <vector>
#include "generalPhysicsAndMath.h"
#include "NJLDimensionfulCouplings.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonPropagators.h"
using namespace std;



class SU3NJL3DCutoffVacuum
{
private:
	SU3NJL3DCutoffParameters parametersNJL;

	double upQuarkEffectiveMass = 0.0/0.0;
	double downQuarkEffectiveMass = 0.0/0.0;
	double strangeQuarkEffectiveMass = 0.0/0.0;

public:
	SU3NJL3DCutoffVacuum(){};
	SU3NJL3DCutoffVacuum(void* );
	SU3NJL3DCutoffVacuum(SU3NJL3DCutoffParameters );

	SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };

	double getUpQuarkEffectiveMass(){ return upQuarkEffectiveMass; };
	double getDownQuarkEffectiveMass(){ return downQuarkEffectiveMass; };
	double getStrangeQuarkEffectiveMass(){ return strangeQuarkEffectiveMass; };

	//gap equations
	void solve(double , MultiRootFindingMethod , double , double , double );
	bool testSolution(double );

	//thermodynamics
	double calculatePressure();
	double calculateEnergyDensity();
	double calculateEntropyDensity(){ return 0.0; }
	double calculateVacuumPressureElectrons(double );

	//meson properties
	SU3NJL3DCutoffMeson calculateMesonMassAndWidth(mesonState , double , MultiRootFindingMethod , double , double );

private:
	void setUpQuarkEffectiveMass(double upQuarkEffectiveMassAux){ upQuarkEffectiveMass = upQuarkEffectiveMassAux; };
	void setDownQuarkEffectiveMass(double downQuarkEffectiveMassAux){ downQuarkEffectiveMass = downQuarkEffectiveMassAux; };
	void setStrangeQuarkEffectiveMass(double strangeQuarkEffectiveMassAux){ strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux; };
};


int SU3NJLGapEquationsVacuum(const gsl_vector *, void *, gsl_vector *);


#endif
#ifndef SU3NJL3DCUTOFFFIXEDCHEMPOTTEMP_H
#define SU3NJL3DCUTOFFFIXEDCHEMPOTTEMP_H

#include <vector>
#include "physics_utils/distribution_functions.h"
#include "njl_model/NJLDimensionfulCouplings.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"

using namespace std;


double linearFit(double , double , double , double , double );


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
	mesonState mesonID;

public:
	SU3NJL3DCutoffFixedChemPotTemp(){};
	SU3NJL3DCutoffFixedChemPotTemp(void* );
	SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffParameters , double , double , double , double );
	SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffParameters , double , double , double );

	SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };

	double getTemperature(){ return temperature; };
	void setTemperature(double temperatureAux){ temperature = temperatureAux; };

	double getUpQuarkEffectiveMass(){ return upQuarkEffectiveMass; };
	double getDownQuarkEffectiveMass(){ return downQuarkEffectiveMass; };
	double getStrangeQuarkEffectiveMass(){ return strangeQuarkEffectiveMass; };
	double getUpQuarkChemicalPotential(){ return upQuarkChemicalPotential; };
	double getDownQuarkChemicalPotential(){ return downQuarkChemicalPotential; };
	double getStrangeQuarkChemicalPotential(){ return strangeQuarkChemicalPotential; };

	mesonState getMesonID(){ return mesonID; };
	void setMesonID(mesonState mesonIDAux){ mesonID = mesonIDAux; };

	//gap equations
	void solve(double , MultiRootFindingMethod , double , double , double );
	bool testSolution(double );

	//meson properties
    SU3NJL3DCutoffMeson calculateMesonMassAndWidth(mesonState , double , MultiRootFindingMethod , double , double );
    void findNondiagonalMesonMottTemperature(mesonState , double , MultiRootFindingMethod , double , double , double , double );

	static void evaluateCrossSectionsEqualLightMasses(
		SU3NJL3DCutoffParameters& , 
		double , 
		MultiRootFindingMethod , 
		double , 
		double , 
		double ,
		int ,
		double ,
		MultiRootFindingMethod ,
		double ,
		int ,
		double ,
		MultiRootFindingMethod,
		double ,
		bool ,
		double ,
		int ,
		int
	);
	
private:
	void setUpQuarkEffectiveMass(double upQuarkEffectiveMassAux){ upQuarkEffectiveMass = upQuarkEffectiveMassAux; };
	void setDownQuarkEffectiveMass(double downQuarkEffectiveMassAux){ downQuarkEffectiveMass = downQuarkEffectiveMassAux; };
	void setStrangeQuarkEffectiveMass(double strangeQuarkEffectiveMassAux){ strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux; };
	void setUpQuarkChemicalPotential(double upQuarkChemicalPotentialAux){ upQuarkChemicalPotential = upQuarkChemicalPotentialAux; };
	void setDownQuarkChemicalPotential(double downQuarkChemicalPotentialAux){ downQuarkChemicalPotential = downQuarkChemicalPotentialAux; };
	void setStrangeQuarkChemicalPotential(double strangeQuarkChemicalPotentialAux){ strangeQuarkChemicalPotential = strangeQuarkChemicalPotentialAux; };
};


int SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature(const gsl_vector *, void *, gsl_vector *);

vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(SU3NJL3DCutoffVacuum , double , int , double , MultiRootFindingMethod );

vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromLowToHighTemperatureAtZeroChemicalPotential(SU3NJL3DCutoffFixedChemPotTemp , double , int , double , MultiRootFindingMethod );

std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromLowToHighTemperature(SU3NJL3DCutoffFixedChemPotTemp , double , int , double , MultiRootFindingMethod );

vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromFiniteTemperatureToFiniteChemicalPotential(SU3NJL3DCutoffFixedChemPotTemp , double , int , double , MultiRootFindingMethod );

vector<SU3NJL3DCutoffMeson> mesonPropertiesFromVacuumToFiniteTemperatureAtZeroChemicalPotential(SU3NJL3DCutoffVacuum , vector<SU3NJL3DCutoffFixedChemPotTemp> , mesonState , double , MultiRootFindingMethod , double , double );

int SU3NJL3DCutoffNondiagonalMesonMottTemperatureFixedChemicalPotentials(const gsl_vector *, void *, gsl_vector *);

SU3NJL3DCutoffFixedChemPotTemp nondiagonalMesonMeltingPoint(SU3NJL3DCutoffVacuum , vector<SU3NJL3DCutoffFixedChemPotTemp> , mesonState , double , MultiRootFindingMethod , double , double );

void evaluateCrossSectionsPaperWithKlevanskyParameterSet(double , double , int , int);

void someVacuumAndThermalPropertiesKlevanskyParameterSet();


#endif
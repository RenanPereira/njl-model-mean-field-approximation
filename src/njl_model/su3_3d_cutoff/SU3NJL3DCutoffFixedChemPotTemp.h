#ifndef SU3NJL3DCUTOFFFIXEDCHEMPOTTEMP_H
#define SU3NJL3DCUTOFFFIXEDCHEMPOTTEMP_H

#include "physics_utils/distribution_functions.h"
#include "njl_model/NJLDimensionfulCouplings.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"

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

	double pressure = 0.0/0.0;
	double energyDensity = 0.0/0.0;
	double entropyDensity = 0.0/0.0;

public:
	SU3NJL3DCutoffFixedChemPotTemp(){}
	SU3NJL3DCutoffFixedChemPotTemp(void* );
	SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffParameters , double , double , double , double );
	SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffParameters , double , double , double );
	SU3NJL3DCutoffFixedChemPotTemp(SU3NJL3DCutoffVacuum &);

	SU3NJL3DCutoffParameters getParametersNJL() const { return parametersNJL; }

	double getTemperature() const { return temperature; }
	void setTemperature(double temperatureAux){ temperature = temperatureAux; }

	double getUpQuarkEffectiveMass() const { return upQuarkEffectiveMass; }
	double getDownQuarkEffectiveMass() const { return downQuarkEffectiveMass; }
	double getStrangeQuarkEffectiveMass() const { return strangeQuarkEffectiveMass; }
	double getUpQuarkChemicalPotential() const { return upQuarkChemicalPotential; }
	double getDownQuarkChemicalPotential() const { return downQuarkChemicalPotential; }
	double getStrangeQuarkChemicalPotential() const { return strangeQuarkChemicalPotential; }

	mesonState getMesonID() const { return mesonID; }
	void setMesonID(mesonState mesonIDAux){ mesonID = mesonIDAux; }

	//gap equations
	void solve(double , MultiRootFindingMethod , double , double , double );
	bool testSolution(double );

	//meson properties
    SU3NJL3DCutoffMeson calculateMesonMassAndWidth(
		mesonState , 
		double , 
		MultiRootFindingMethod , 
		double , 
		double 
	);
    void findNondiagonalMesonMottTemperature(
		mesonState , 
		double , 
		MultiRootFindingMethod , 
		double , 
		double , 
		double , 
		double 
	);

	static void evaluateIsospinSymmetricCrossSections(
		SU3NJL3DCutoffParameters& , 
		double , 
		MultiRootFindingMethod , 
		double , 
		double , 
		double , 
		double , 
		int , 
		double , 
		MultiRootFindingMethod , 
		double , 
		int , 
		double , 
		MultiRootFindingMethod , 
		double , 
		bool , 
		double , 
		int , 
		int 
	);

	double calculatePressure(double );
	double calculateEnergyDensity(double );
	double calculateEntropyDensity();

	void setPressure(double pressureVacuum){ pressure = calculatePressure(pressureVacuum); }
	void setEnergyDensity(double energyVacuum){ energyDensity = calculateEnergyDensity(energyVacuum); }
	void setEntropyDensity(){ entropyDensity = calculateEntropyDensity(); }

	double getPressure() const { return pressure; }
	double getEnergyDensity() const { return energyDensity; }
	double getEntropyDensity() const { return entropyDensity; }

	static void evaluateInMediumMassesAndThermodynamics(
		SU3NJL3DCutoffParameters& , 
		double , 
		MultiRootFindingMethod , 
		double , 
		double , 
		double , 
		double ,
		double , 
		int , 
		double , 
		MultiRootFindingMethod 
	);

	static std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromFiniteTemperatureToFiniteChemicalPotential(
		SU3NJL3DCutoffFixedChemPotTemp , 
		double , 
		int , 
		double , 
		MultiRootFindingMethod 
	);

	static std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromLowToHighTemperatureAtZeroChemicalPotential(
		SU3NJL3DCutoffFixedChemPotTemp , 
		double , 
		int , 
		double , 
		MultiRootFindingMethod 
	);

	static std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
		SU3NJL3DCutoffVacuum , 
		double ,
		double , 
		int , 
		double , 
		MultiRootFindingMethod 
	);

	static std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveFromLowToHighTemperature(
		SU3NJL3DCutoffFixedChemPotTemp , 
		double , 
		int , 
		double , 
		MultiRootFindingMethod 
	);

	static std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveToTemperature(
		SU3NJL3DCutoffFixedChemPotTemp , 
		double , 
		int , 
		double , 
		MultiRootFindingMethod 
	);

	static vector<SU3NJL3DCutoffFixedChemPotTemp> solveToChemicalPotentialSymmetric(
		SU3NJL3DCutoffFixedChemPotTemp , 
		double , 
		int , 
		double , 
		MultiRootFindingMethod 
	);

	static std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveVacuumToChemicalPotential(
		SU3NJL3DCutoffVacuum& , 
		double , 
		int , 
		double , 
		MultiRootFindingMethod 
	);
	
	static std::vector<SU3NJL3DCutoffFixedChemPotTemp> solveVacuumToTemperature(
		SU3NJL3DCutoffVacuum& , 
		double , 
		int , 
		double , 
		MultiRootFindingMethod 
	);

	static void computeThermoFixedChemPotTrajectory(
		SU3NJL3DCutoffParameters& , 
		double ,
		MultiRootFindingMethod ,
		double ,
		double ,
		double ,
		double ,
		int ,
		double ,
		MultiRootFindingMethod ,
		double ,
		int ,
		double ,
		MultiRootFindingMethod ,
		std::string& 
	);

	static void computeThermoFixedTemperatureTrajectory(
		SU3NJL3DCutoffParameters& , 
		double ,
		MultiRootFindingMethod ,
		double ,
		double ,
		double ,
		double ,
		int ,
		double ,
		MultiRootFindingMethod ,
		double ,
		int ,
		double ,
		MultiRootFindingMethod ,
		string& 
	);
	
private:
	void setUpQuarkEffectiveMass(double upQuarkEffectiveMassAux){ upQuarkEffectiveMass = upQuarkEffectiveMassAux; }
	void setDownQuarkEffectiveMass(double downQuarkEffectiveMassAux){ downQuarkEffectiveMass = downQuarkEffectiveMassAux; }
	void setStrangeQuarkEffectiveMass(double strangeQuarkEffectiveMassAux){ strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux; }
	
	void setUpQuarkChemicalPotential(double upQuarkChemicalPotentialAux){ upQuarkChemicalPotential = upQuarkChemicalPotentialAux; }
	void setDownQuarkChemicalPotential(double downQuarkChemicalPotentialAux){ downQuarkChemicalPotential = downQuarkChemicalPotentialAux; }
	void setStrangeQuarkChemicalPotential(double strangeQuarkChemicalPotentialAux){ strangeQuarkChemicalPotential = strangeQuarkChemicalPotentialAux; }
};

int SU3NJL3DCutoffGapEquationsFixedChemicalPotentialsTemperature(const gsl_vector *, void *, gsl_vector *);

std::vector<SU3NJL3DCutoffMeson> mesonPropertiesFromVacuumToFiniteTemperatureAtZeroChemicalPotential(
	SU3NJL3DCutoffVacuum , 
	std::vector<SU3NJL3DCutoffFixedChemPotTemp> , 
	mesonState , 
	double , 
	MultiRootFindingMethod , 
	double , 
	double 
);

int SU3NJL3DCutoffNondiagonalMesonMottTemperatureFixedChemicalPotentials(const gsl_vector *, void *, gsl_vector *);

SU3NJL3DCutoffFixedChemPotTemp nondiagonalMesonMeltingPoint(
	SU3NJL3DCutoffVacuum , 
	std::vector<SU3NJL3DCutoffFixedChemPotTemp> , 
	mesonState , 
	double , 
	MultiRootFindingMethod , 
	double , 
	double 
);

void writeSolutionsToFile(std::vector<SU3NJL3DCutoffFixedChemPotTemp> , std::string );

void calculateThermodynamics(SU3NJL3DCutoffVacuum& , std::vector<SU3NJL3DCutoffFixedChemPotTemp>& );

#endif

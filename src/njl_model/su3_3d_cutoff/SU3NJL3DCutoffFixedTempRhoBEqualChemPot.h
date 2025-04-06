#ifndef SU3NJL3DCUTOFFEQUALCHEMPOTFIXEDTEMPRHOB_H
#define SU3NJL3DCUTOFFEQUALCHEMPOTFIXEDTEMPRHOB_H

#include <vector>
#include "gsl_wrapper/root_solver_gsl.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h"


class SU3NJL3DCutoffFixedTempRhoBEqualChemPot
{
public:
	struct ChiralTransitionPoint {
		double mU_broken, mD_broken, mS_broken;
		double effCP_broken;
		double mU_restored, mD_restored, mS_restored;
		double effCP_restored;
	};


private:
	SU3NJL3DCutoffParameters parametersNJL;

	double temperature = 0.0;
	double baryonDensity = 0.0;

	double upQuarkEffectiveMass = 0.0/0.0;
	double downQuarkEffectiveMass = 0.0/0.0;
	double strangeQuarkEffectiveMass = 0.0/0.0;
	
	double quarkEffectiveChemicalPotential = 0.0/0.0;

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

	//Thermodynamics
	double pressure = 0.0/0.0;
	double energyDensity = 0.0/0.0;
	double entropyDensity = 0.0/0.0;


public:
	SU3NJL3DCutoffFixedTempRhoBEqualChemPot(){};
	SU3NJL3DCutoffFixedTempRhoBEqualChemPot(void* );
	SU3NJL3DCutoffFixedTempRhoBEqualChemPot(SU3NJL3DCutoffParameters );
	SU3NJL3DCutoffFixedTempRhoBEqualChemPot(SU3NJL3DCutoffParameters , double );
	SU3NJL3DCutoffFixedTempRhoBEqualChemPot(SU3NJL3DCutoffParameters , double , double );
	SU3NJL3DCutoffFixedTempRhoBEqualChemPot(SU3NJL3DCutoffParameters , double , double , double , double , double , double );
	SU3NJL3DCutoffFixedTempRhoBEqualChemPot(SU3NJL3DCutoffVacuum );

	SU3NJL3DCutoffParameters getParametersNJL(){ return parametersNJL; };

	double getTemperature(){ return temperature; };
	double getBaryonDensity(){ return baryonDensity; };

	void setPressure(double pressureVacuum){ pressure = calculatePressure(pressureVacuum); };
	void setEnergyDensity(double energyVacuum){ energyDensity = calculateEnergyDensity(energyVacuum); };
	void setEntropyDensity(){ entropyDensity = calculateEntropyDensity(); };

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

	//baryon chemical potential
	void setBaryonChemicalPotential(){ baryonChemicalPotential = getUpQuarkChemicalPotential() + 2.0*getDownQuarkChemicalPotential(); };
	double getBaryonChemicalPotential(){ return baryonChemicalPotential; }
	
	void setSigmasDensitiesChemicalPotentials(double , double , double , double , double , double );
	void setSigmasDensitiesChemicalPotentials();
	
	void setUpQuarkEffectiveMass(double upQuarkEffectiveMassAux){ upQuarkEffectiveMass = upQuarkEffectiveMassAux; };
	void setDownQuarkEffectiveMass(double downQuarkEffectiveMassAux){ downQuarkEffectiveMass = downQuarkEffectiveMassAux; };
	void setStrangeQuarkEffectiveMass(double strangeQuarkEffectiveMassAux){ strangeQuarkEffectiveMass = strangeQuarkEffectiveMassAux; };
	void setQuarkEffectiveChemicalPotential(double quarkEffectiveChemicalPotentialAux){ quarkEffectiveChemicalPotential = quarkEffectiveChemicalPotentialAux; };

	static ChiralTransitionPoint calculateChiralTransitionPoint(
		double ,
		const ChiralTransitionPoint& ,
		const SU3NJL3DCutoffParameters& ,
		double ,
		MultiRootFindingMethod );
		static vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot::ChiralTransitionPoint> calculateFirstOrderLine(vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> , double , MultiRootFindingMethod );


private:
	void setTemperature(double temperatureAux){ temperature = temperatureAux; };
	void setBaryonDensity(double baryonDensityAux){ baryonDensity = baryonDensityAux; };
	
	void setUpQuarkSigma(double upQuarkSigmaAux){ upQuarkSigma = upQuarkSigmaAux; };
	void setDownQuarkSigma(double downQuarkSigmaAux){ downQuarkSigma = downQuarkSigmaAux; };
	void setStrangeQuarkSigma(double strangeQuarkSigmaAux){ strangeQuarkSigma = strangeQuarkSigmaAux; };

	void setUpQuarkDensity(double upQuarkDensityAux){ upQuarkDensity = upQuarkDensityAux; };
	void setDownQuarkDensity(double downQuarkDensityAux){ downQuarkDensity = downQuarkDensityAux; };
	void setStrangeQuarkDensity(double strangeQuarkDensityAux){ strangeQuarkDensity = strangeQuarkDensityAux; };

	void setUpQuarkChemicalPotential(double upQuarkChemicalPotentialAux){ upQuarkChemicalPotential = upQuarkChemicalPotentialAux; };
	void setDownQuarkChemicalPotential(double downQuarkChemicalPotentialAux){ downQuarkChemicalPotential = downQuarkChemicalPotentialAux; };
	void setStrangeQuarkChemicalPotential(double strangeQuarkChemicalPotentialAux){ strangeQuarkChemicalPotential = strangeQuarkChemicalPotentialAux; };
};


int SU3NJL3DCutoffGapEquationsEqualChemicalPotentialFixedTemperatureBaryonDensity(const gsl_vector *, void *, gsl_vector *);

int SU3NJL3DCutoffEqualChemPotFixedTempChiralTransitionPoint(const gsl_vector *, void *, gsl_vector *);

vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> 
solveFromVacuumToFiniteBaryonDensity(SU3NJL3DCutoffVacuum , 
                                     double , double , int , 
                                     double , MultiRootFindingMethod , bool);

void writeSolutionsToFile(vector<SU3NJL3DCutoffFixedTempRhoBEqualChemPot> , string , bool );


#endif
SHELL := /bin/bash
#CXX = g++ -std=c++11 -fopenmp
#CXX = g++ -O3 -std=c++11 -fopenmp
CXX = g++ -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp


DEPS = src/OneVariableFunction.h src/Integration1DimNewtonCotes.h src/UnitaryGroup3Dimensions.h \
       src/rootSolverGSL.h src/Interpolation1DimGSL.h src/Integration1DimGSL.h src/ComplexSquareMatrixGSL.h \
       src/generalPhysicsAndMath.h src/OneFermionLineIntegral.h src/TwoFermionLineIntegral.h \
       src/NJLDimensionfulCouplings.h \
       src/SU3NJL3DCutoff.h \
       src/SU3NJL3DCutoffVacuum.h \
       src/SU3NJL3DCutoffFixedChemPotTemp.h \
       src/SU3NJL3DCutoffEqualChemPotFixedTempRhoB.h \
       src/SU3NJL3DCutoffBetaEqFixedTempRhoB.h \
       src/SU3NJL3DCutoffMesonProjectors.h src/SU3NJL3DCutoffMesonPropagators.h \
       src/SU3NJL3DCutoffDifferentialCrossSections.h src/SU3NJL3DCutoffCrossSections.h src/SU3NJL3DCutoffIntegratedCrossSections.h 


OBJ = obj/main.o \
      obj/OneVariableFunction.o  obj/Integration1DimNewtonCotes.o obj/UnitaryGroup3Dimensions.o \
      obj/rootSolverGSL.o obj/Interpolation1DimGSL.o obj/Integration1DimGSL.o obj/ComplexSquareMatrixGSL.o \
      obj/generalPhysicsAndMath.o obj/OneFermionLineIntegral.o obj/TwoFermionLineIntegral.o \
      obj/NJLDimensionfulCouplings.o \
      obj/SU3NJL3DCutoff.o \
      obj/SU3NJL3DCutoffVacuum.o \
      obj/SU3NJL3DCutoffFixedChemPotTemp.o \
      obj/SU3NJL3DCutoffEqualChemPotFixedTempRhoB.o \
      obj/SU3NJL3DCutoffBetaEqFixedTempRhoB.o \
      obj/SU3NJL3DCutoffMesonProjectors.o obj/SU3NJL3DCutoffMesonPropagators.o \
      obj/SU3NJL3DCutoffDifferentialCrossSections.o obj/SU3NJL3DCutoffCrossSections.o obj/SU3NJL3DCutoffIntegratedCrossSections.o


obj/%.o: src/%.cpp $(DEPS)
	$(CXX) $< -Iinclude -c -o $@

run: $(OBJ) 
	$(CXX) $(OBJ) -Iinclude -L/usr/local/lib -lgsl -lgslcblas -o bin/nambuJonaLasinioModel.out

clean: 
	rm -f bin/nambuJonaLasinioModel.out $(OBJ)



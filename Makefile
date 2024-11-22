SHELL := /bin/bash
#CXX = g++ -std=c++11 -fopenmp
#CXX = g++ -O3 -std=c++11 -fopenmp
CXX = g++ -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp
#CXX = g++ -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp

INCLUDE_DIRS = -Isrc -Isrc/ini_file_parser -Isrc/math_utils -Isrc/gsl_wrapper

DEPS = src/math_utils/OneVariableFunction.h \
       src/Integration1DimNewtonCotes.h \
       src/UnitaryGroup3Dimensions.h \
       src/gsl_wrapper/root_solver_gsl.h \
       src/gsl_wrapper/Interpolation1DimGSL.h \
       src/gsl_wrapper/Integration1DimGSL.h \
       src/gsl_wrapper/ComplexSquareMatrixGSL.h \
       src/generalPhysicsAndMath.h \
       src/ini_file_parser/IniFileParser.h \
       src/OneFermionLineIntegral.h \
       src/TwoFermionLineIntegral.h \
       src/NJLDimensionfulCouplings.h \
       src/SU3NJL3DCutoff.h \
       src/SU3NJL3DCutoffVacuum.h \
       src/SU3NJL3DCutoffFixedChemPotTemp.h \
       src/SU3NJL3DCutoffEqualChemPotFixedTempRhoB.h \
       src/SU3NJL3DCutoffBetaEqFixedTempRhoB.h \
       src/SU3NJL3DCutoffMesonProjectors.h \
       src/SU3NJL3DCutoffMesonPropagators.h \
       src/SU3NJL3DCutoffDifferentialCrossSections.h \
       src/SU3NJL3DCutoffCrossSections.h \
       src/SU3NJL3DCutoffIntegratedCrossSections.h 


OBJ = obj/main.o \
      obj/math_utils/OneVariableFunction.o \
      obj/Integration1DimNewtonCotes.o \
      obj/UnitaryGroup3Dimensions.o \
      obj/gsl_wrapper/root_solver_gsl.o \
      obj/gsl_wrapper/Interpolation1DimGSL.o \
      obj/gsl_wrapper/Integration1DimGSL.o \
      obj/gsl_wrapper/ComplexSquareMatrixGSL.o \
      obj/generalPhysicsAndMath.o \
      obj/ini_file_parser/IniFileParser.o \
      obj/OneFermionLineIntegral.o \
      obj/TwoFermionLineIntegral.o \
      obj/NJLDimensionfulCouplings.o \
      obj/SU3NJL3DCutoff.o \
      obj/SU3NJL3DCutoffVacuum.o \
      obj/SU3NJL3DCutoffFixedChemPotTemp.o \
      obj/SU3NJL3DCutoffEqualChemPotFixedTempRhoB.o \
      obj/SU3NJL3DCutoffBetaEqFixedTempRhoB.o \
      obj/SU3NJL3DCutoffMesonProjectors.o \
      obj/SU3NJL3DCutoffMesonPropagators.o \
      obj/SU3NJL3DCutoffDifferentialCrossSections.o \
      obj/SU3NJL3DCutoffCrossSections.o \
      obj/SU3NJL3DCutoffIntegratedCrossSections.o


obj/%.o: src/%.cpp $(DEPS)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/ini_file_parser/%.o: src/ini_file_parser/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/math_utils/%.o: src/math_utils/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/gsl_wrapper/%.o: src/gsl_wrapper/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

run: $(OBJ) 
	$(CXX) $(OBJ) $(INCLUDE_DIRS) -L/usr/local/lib -lgsl -lgslcblas -o bin/nambuJonaLasinioModel.out

clean: 
	rm -f bin/nambuJonaLasinioModel.out $(OBJ)



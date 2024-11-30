SHELL := /bin/bash
#CXX = g++ -std=c++11 -fopenmp
#CXX = g++ -O3 -std=c++11 -fopenmp
CXX = g++ -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp
#CXX = g++ -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp

INCLUDE_DIRS = -Isrc

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
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoff.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffEqualChemPotFixedTempRhoB.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffBetaEqFixedTempRhoB.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonProjectors.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonPropagators.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffDifferentialCrossSections.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffCrossSections.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.h 


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
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoff.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffEqualChemPotFixedTempRhoB.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffBetaEqFixedTempRhoB.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonProjectors.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonPropagators.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffDifferentialCrossSections.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffCrossSections.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.o


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

obj/njl_model/su3_3d_cutoff/%.o: src/njl_model/su3_3d_cutoff/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@


run: $(OBJ) 
	$(CXX) $(OBJ) $(INCLUDE_DIRS) -L/usr/local/lib -lgsl -lgslcblas -o bin/nambuJonaLasinioModel.out

clean: 
	rm -f bin/nambuJonaLasinioModel.out $(OBJ)



SHELL := /bin/bash
#CXX = g++ -std=c++11 -fopenmp
#CXX = g++ -O3 -std=c++11 -fopenmp
CXX = g++ -O3 -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp
#CXX = g++ -Wall -Wextra -Wfloat-equal -Wundef -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -std=c++11 -fopenmp


INCLUDE_DIRS = -Isrc


DEPS = src/command_line_processor.h \
       src/math_utils/OneVariableFunction.h \
       src/math_utils/useful_functions.h \
       src/integration_methods/Integration1DimNewtonCotes.h \
       src/group_theory/UnitaryGroup3Dimensions.h \
       src/gsl_wrapper/root_solver_gsl.h \
       src/gsl_wrapper/InterpolationGSL1Dim.h \
       src/gsl_wrapper/Integration1DimGSL.h \
       src/gsl_wrapper/ComplexSquareMatrixGSL.h \
       src/utils/format_utils.h \
       src/physics_utils/distribution_functions.h \
       src/physics_utils/physical_constants.h \
       src/ini_file_parser/IniFileParser.h \
       src/njl_model/NJLDimensionlessCouplings.h \
       src/njl_model/NJLDimensionfulCouplings.h \
       src/njl_model/njl_regularization_schemes.h \
       src/njl_model/n_fermion_line_integrals/one_fermion_line_integral_3d_cutoff.h \
       src/njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff.h \
       src/njl_model/n_fermion_line_integrals/KlevanskyB0Integral3DCutoffFileParser.h \
       src/njl_model/n_fermion_line_integrals/n_fermion_line_integrals_calculator.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoff.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedTempRhoBEqualChemPot.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffBetaEqFixedTempRhoB.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonProjectors.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonPropagators.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffDifferentialCrossSections.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffCrossSections.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.h \
       src/njl_model/su3_3d_cutoff/SU3NJL3DCutoffCalculator.h


OBJ = obj/main.o \
      obj/command_line_processor.o \
      obj/math_utils/OneVariableFunction.o \
      obj/math_utils/useful_functions.o \
      obj/integration_methods/Integration1DimNewtonCotes.o \
      obj/group_theory/UnitaryGroup3Dimensions.o \
      obj/gsl_wrapper/root_solver_gsl.o \
      obj/gsl_wrapper/InterpolationGSL1Dim.o \
      obj/gsl_wrapper/Integration1DimGSL.o \
      obj/gsl_wrapper/ComplexSquareMatrixGSL.o \
      obj/utils/format_utils.o \
      obj/physics_utils/distribution_functions.o \
      obj/ini_file_parser/IniFileParser.o \
      obj/njl_model/NJLDimensionfulCouplings.o \
      obj/njl_model/njl_regularization_schemes.o \
      obj/njl_model/n_fermion_line_integrals/one_fermion_line_integral_3d_cutoff.o \
      obj/njl_model/n_fermion_line_integrals/two_fermion_line_integral_3d_cutoff.o \
      obj/njl_model/n_fermion_line_integrals/KlevanskyB0Integral3DCutoffFileParser.o \
      obj/njl_model/n_fermion_line_integrals/n_fermion_line_integrals_calculator.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoff.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffVacuum.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedTempRhoBEqualChemPot.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffBetaEqFixedTempRhoB.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonProjectors.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonPropagators.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffDifferentialCrossSections.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffCrossSections.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffIntegratedCrossSections.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffFileParser.o \
      obj/njl_model/su3_3d_cutoff/SU3NJL3DCutoffCalculator.o


obj/%.o: src/%.cpp $(DEPS)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/math_utils/%.o: src/math_utils/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/integration_methods/%.o: src/integration_methods/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/group_theory/%.o: src/group_theory/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/gsl_wrapper/%.o: src/gsl_wrapper/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/utils/%.o: src/utils/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/physics_utils/%.o: src/physics_utils/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/ini_file_parser/%.o: src/ini_file_parser/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/njl_model/n_fermion_line_integrals/%.o: src/njl_model/n_fermion_line_integrals/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@

obj/njl_model/su3_3d_cutoff/%.o: src/njl_model/su3_3d_cutoff/%.cpp $(DEPS)
	@mkdir -p $(dir $@)
	$(CXX) $(INCLUDE_DIRS) -c $< -o $@


run: $(OBJ) 
	$(CXX) $(OBJ) $(INCLUDE_DIRS) -L/usr/local/lib -lgsl -lgslcblas -o bin/nambuJonaLasinioModel.out

clean: 
	rm -f bin/nambuJonaLasinioModel.out $(OBJ)



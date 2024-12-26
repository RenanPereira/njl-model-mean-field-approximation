#ifndef SU3NJL3DCUTOFFCALCULATOR_H
#define SU3NJL3DCUTOFFCALCULATOR_H

#include "ini_file_parser/IniFileParser.h"
#include "gsl_wrapper/root_solver_gsl.h"
#include "SU3NJL3DCutoff.h"


void evaluateSU3NJL3DCutoffVacuumMasses(SU3NJL3DCutoffParameters& ,
                                        double ,
                                        MultiRootFindingMethod ,
                                        double , double , double );

void evaluateSU3NJL3DCutoffVacuumMasses(const IniFileParser& );

#endif

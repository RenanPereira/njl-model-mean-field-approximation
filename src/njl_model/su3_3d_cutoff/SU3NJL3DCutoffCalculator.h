#ifndef SU3NJL3DCUTOFFCALCULATOR_H
#define SU3NJL3DCUTOFFCALCULATOR_H

#include "ini_file_parser/IniFileParser.h"
#include "gsl_wrapper/root_solver_gsl.h"
#include "SU3NJL3DCutoff.h"

namespace SU3NJL3DCutoffCalculator
{
    void evaluateVacuumMasses(
        SU3NJL3DCutoffParameters& ,                                
        double ,                                
        MultiRootFindingMethod ,                                
        double , double , double 
    );
    void evaluateVacuumMasses(const IniFileParser& );

    void evaluateFirstOrderLine(
        SU3NJL3DCutoffParameters& ,                                
        double ,                                
        MultiRootFindingMethod ,                                
        double , 
        double , 
        double ,
        double , 
        double, 
        int ,
        double ,
        MultiRootFindingMethod ,
        bool ,
        double , 
        MultiRootFindingMethod ,
        double ,
        double 
    );
}

#endif

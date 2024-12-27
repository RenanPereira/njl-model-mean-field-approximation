#include "math_utils/useful_functions.h"

//Function that returns the sign of x
double sign(double x)
{   
    if ( x>0 ){ return +1.0; }
    else
    {
        if ( x<0 ){ return -1.0; }
        else{ return 0.0; }
    }
}

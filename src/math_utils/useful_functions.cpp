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

//linear interpolation and extrapolation
double linearFit(double x, double x1, double y1, double x2, double y2)
{
    double m = ( y1 - y2 )/( x1 - x2 );
    double b = 0.5*( ( y1 + y2 ) - m*( x1 + x2 ) );
    double y = m*x + b;

    return y;
}

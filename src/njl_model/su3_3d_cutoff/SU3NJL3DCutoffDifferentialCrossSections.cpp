#include <cmath>
#include <iostream>
#include <gsl/gsl_complex_math.h>
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonPropagators.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffDifferentialCrossSections.h"



//In this file we have functions to calculate quark-quark and quark-antiquark differential cross sections
//quark-quark: 
//UD->UD, DU->DU, US->US, SU->SU, DS->DS, SD->SD, UU->UU, DD->DD, SS->SS
//quark-antiquark: 
//UDBar->UDBar, USBar->USBar, DSBar->DSBar,
//DUBar->DUBar, SUBar->SUBar, SDBar->SDBar,
//UUBar->UUBar, UUBar->DDBar, UUBar->SSBar,
//DDBar->UUBar, DDBar->DDBar, DDBar->SSBar,
//SSBar->UUBar, SSBar->DDBar, SSBar->SSBar
//antiquark-antiquark: 
//UBarUBar->UBarUBar, UBarDBar->UBarDBar, UBarSBar->UBarSBar, SBarSBar->SBarSBar


string toString(scatteringProcess process)
{ 
    if      ( process==UDUD ){ return "UDUD"; }
    else if ( process==DUDU ){ return "DUDU"; }
    else if ( process==USUS ){ return "USUS"; }
    else if ( process==SUSU ){ return "SUSU"; }
    else if ( process==DSDS ){ return "DSDS"; }
    else if ( process==SDSD ){ return "SDSD"; }
    else if ( process==UUUU ){ return "UUUU"; }
    else if ( process==DDDD ){ return "DDDD"; }
    else if ( process==SSSS ){ return "SSSS"; }
    else if ( process==UDBarUDBar ){ return "UDBarUDBar"; }
    else if ( process==USBarUSBar ){ return "USBarUSBar"; }
    else if ( process==DSBarDSBar ){ return "DSBarDSBar"; }
    else if ( process==DUBarDUBar ){ return "DUBarDUBar"; }
    else if ( process==SUBarSUBar ){ return "SUBarSUBar"; }
    else if ( process==SDBarSDBar ){ return "SDBarSDBar"; }
    else if ( process==UUBarUUBar ){ return "UUBarUUBar"; }
    else if ( process==UUBarDDBar ){ return "UUBarDDBar"; }
    else if ( process==UUBarSSBar ){ return "UUBarSSBar"; }
    else if ( process==DDBarUUBar ){ return "DDBarUUBar"; }
    else if ( process==DDBarDDBar ){ return "DDBarDDBar"; }
    else if ( process==DDBarSSBar ){ return "DDBarSSBar"; }
    else if ( process==SSBarUUBar ){ return "SSBarUUBar"; }
    else if ( process==SSBarDDBar ){ return "SSBarDDBar"; }
    else if ( process==SSBarSSBar ){ return "SSBarSSBar"; }
    else if ( process==UBarUBarUBarUBar ){ return "UBarUBarUBarUBar"; }
    else if ( process==DBarDBarDBarDBar ){ return "DBarDBarDBarDBar"; }
    else if ( process==SBarSBarSBarSBar ){ return "SBarSBarSBarSBar"; }
    else if ( process==UBarDBarUBarDBar ){ return "UBarDBarUBarDBar"; }
    else if ( process==UBarSBarUBarSBar ){ return "UBarSBarUBarSBar"; }

    return "0";
}


//function that given the process selects the incoming and outgoing masses
void inOutMassesGivenScatteringProcess(double effMU, double effMD, double effMS, scatteringProcess process, double& m1, double& m2, double& m3, double& m4)
{   
    if      ( process==UDUD ){ m1 = effMU; m2 = effMD; m3 = effMU; m4 = effMD; }
    else if ( process==DUDU ){ m1 = effMD; m2 = effMU; m3 = effMD; m4 = effMU; }
    else if ( process==USUS ){ m1 = effMU; m2 = effMS; m3 = effMU; m4 = effMS; }
    else if ( process==SUSU ){ m1 = effMS; m2 = effMU; m3 = effMS; m4 = effMU; }
    else if ( process==DSDS ){ m1 = effMD; m2 = effMS; m3 = effMD; m4 = effMS; }
    else if ( process==SDSD ){ m1 = effMS; m2 = effMD; m3 = effMS; m4 = effMD; }
    else if ( process==UUUU ){ m1 = effMU; m2 = effMU; m3 = effMU; m4 = effMU; }
    else if ( process==DDDD ){ m1 = effMD; m2 = effMD; m3 = effMD; m4 = effMD; }
    else if ( process==SSSS ){ m1 = effMS; m2 = effMS; m3 = effMS; m4 = effMS; }
    else if ( process==UDBarUDBar ){ m1 = effMU; m2 = effMD; m3 = effMU; m4 = effMD; }
    else if ( process==DUBarDUBar ){ m1 = effMD; m2 = effMU; m3 = effMD; m4 = effMU; }
    else if ( process==USBarUSBar ){ m1 = effMU; m2 = effMS; m3 = effMU; m4 = effMS; }
    else if ( process==SUBarSUBar ){ m1 = effMS; m2 = effMU; m3 = effMS; m4 = effMU; }
    else if ( process==DSBarDSBar ){ m1 = effMD; m2 = effMS; m3 = effMD; m4 = effMS; }
    else if ( process==SDBarSDBar ){ m1 = effMS; m2 = effMD; m3 = effMS; m4 = effMD; }
    else if ( process==UUBarUUBar ){ m1 = effMU; m2 = effMU; m3 = effMU; m4 = effMU; }
    else if ( process==UUBarDDBar ){ m1 = effMU; m2 = effMU; m3 = effMD; m4 = effMD; }
    else if ( process==UUBarSSBar ){ m1 = effMU; m2 = effMU; m3 = effMS; m4 = effMS; }
    else if ( process==DDBarUUBar ){ m1 = effMD; m2 = effMD; m3 = effMU; m4 = effMU; }
    else if ( process==DDBarDDBar ){ m1 = effMD; m2 = effMD; m3 = effMD; m4 = effMD; }
    else if ( process==DDBarSSBar ){ m1 = effMD; m2 = effMD; m3 = effMS; m4 = effMS; }
    else if ( process==SSBarUUBar ){ m1 = effMS; m2 = effMS; m3 = effMU; m4 = effMU; }
    else if ( process==SSBarDDBar ){ m1 = effMS; m2 = effMS; m3 = effMD; m4 = effMD; }
    else if ( process==SSBarSSBar ){ m1 = effMS; m2 = effMS; m3 = effMS; m4 = effMS; }
    else if ( process==UBarUBarUBarUBar ){ m1 = effMU; m2 = effMU; m3 = effMU; m4 = effMU; }
    else if ( process==DBarDBarDBarDBar ){ m1 = effMD; m2 = effMD; m3 = effMD; m4 = effMD; }
    else if ( process==SBarSBarSBarSBar ){ m1 = effMS; m2 = effMS; m3 = effMS; m4 = effMS; }
    else if ( process==UBarDBarUBarDBar ){ m1 = effMU; m2 = effMD; m3 = effMU; m4 = effMD; }
    else if ( process==UBarSBarUBarSBar ){ m1 = effMU; m2 = effMS; m3 = effMU; m4 = effMS; }
}


//function that given the process selects the incoming chemical potentials
void inChemicalPotentialsGivenScatteringProcess(double effCPU, double effCPD, double effCPS, scatteringProcess process, double& cP1, double& cP2)
{   
    if      ( process==UDUD ){ cP1 = effCPU; cP2 = effCPD; }
    else if ( process==DUDU ){ cP1 = effCPD; cP2 = effCPU; }
    else if ( process==USUS ){ cP1 = effCPU; cP2 = effCPS; }
    else if ( process==SUSU ){ cP1 = effCPS; cP2 = effCPU; }
    else if ( process==DSDS ){ cP1 = effCPD; cP2 = effCPS; }
    else if ( process==SDSD ){ cP1 = effCPS; cP2 = effCPD; }
    else if ( process==UUUU ){ cP1 = effCPU; cP2 = effCPU; }
    else if ( process==DDDD ){ cP1 = effCPD; cP2 = effCPD; }
    else if ( process==SSSS ){ cP1 = effCPS; cP2 = effCPS; }
    else if ( process==UDBarUDBar ){ cP1 = effCPU; cP2 = -effCPD; }
    else if ( process==DUBarDUBar ){ cP1 = effCPD; cP2 = -effCPU; }
    else if ( process==USBarUSBar ){ cP1 = effCPU; cP2 = -effCPS; }
    else if ( process==SUBarSUBar ){ cP1 = effCPS; cP2 = -effCPU; }
    else if ( process==DSBarDSBar ){ cP1 = effCPD; cP2 = -effCPS; }
    else if ( process==SDBarSDBar ){ cP1 = effCPS; cP2 = -effCPD; }
    else if ( process==UUBarUUBar ){ cP1 = effCPU; cP2 = -effCPU; }
    else if ( process==UUBarDDBar ){ cP1 = effCPU; cP2 = -effCPU; }
    else if ( process==UUBarSSBar ){ cP1 = effCPU; cP2 = -effCPU; }
    else if ( process==DDBarUUBar ){ cP1 = effCPD; cP2 = -effCPD; }
    else if ( process==DDBarDDBar ){ cP1 = effCPD; cP2 = -effCPD; }
    else if ( process==DDBarSSBar ){ cP1 = effCPD; cP2 = -effCPD; }
    else if ( process==SSBarUUBar ){ cP1 = effCPS; cP2 = -effCPS; }
    else if ( process==SSBarDDBar ){ cP1 = effCPS; cP2 = -effCPS; }
    else if ( process==SSBarSSBar ){ cP1 = effCPS; cP2 = -effCPS; }
    else if ( process==UBarUBarUBarUBar ){ cP1 = -effCPU; cP2 = -effCPU; }
    else if ( process==DBarDBarDBarDBar ){ cP1 = -effCPD; cP2 = -effCPD; }
    else if ( process==SBarSBarSBarSBar ){ cP1 = -effCPS; cP2 = -effCPS; }
    else if ( process==UBarDBarUBarDBar ){ cP1 = -effCPU; cP2 = -effCPD; }
    else if ( process==UBarSBarUBarSBar ){ cP1 = -effCPU; cP2 = -effCPS; }
}


void outChemicalPotentialsGivenScatteringProcess(double effCPU, double effCPD, double effCPS, scatteringProcess process, double& cP3, double& cP4)
{   
    if      ( process==UDUD ){ cP3 = effCPU; cP4 = effCPD; }
    else if ( process==DUDU ){ cP3 = effCPD; cP4 = effCPU; }
    else if ( process==USUS ){ cP3 = effCPU; cP4 = effCPS; }
    else if ( process==SUSU ){ cP3 = effCPS; cP4 = effCPU; }
    else if ( process==DSDS ){ cP3 = effCPD; cP4 = effCPS; }
    else if ( process==SDSD ){ cP3 = effCPS; cP4 = effCPD; }
    else if ( process==UUUU ){ cP3 = effCPU; cP4 = effCPU; }
    else if ( process==DDDD ){ cP3 = effCPD; cP4 = effCPD; }
    else if ( process==SSSS ){ cP3 = effCPS; cP4 = effCPS; }
    else if ( process==UDBarUDBar ){ cP3 = effCPU; cP4 = -effCPD; }
    else if ( process==DUBarDUBar ){ cP3 = effCPD; cP4 = -effCPU; }
    else if ( process==USBarUSBar ){ cP3 = effCPU; cP4 = -effCPS; }
    else if ( process==SUBarSUBar ){ cP3 = effCPS; cP4 = -effCPU; }
    else if ( process==DSBarDSBar ){ cP3 = effCPD; cP4 = -effCPS; }
    else if ( process==SDBarSDBar ){ cP3 = effCPS; cP4 = -effCPD; }
    else if ( process==UUBarUUBar ){ cP3 = effCPU; cP4 = -effCPU; }
    else if ( process==UUBarDDBar ){ cP3 = effCPD; cP4 = -effCPD; }
    else if ( process==UUBarSSBar ){ cP3 = effCPS; cP4 = -effCPS; }
    else if ( process==DDBarUUBar ){ cP3 = effCPU; cP4 = -effCPU; }
    else if ( process==DDBarDDBar ){ cP3 = effCPD; cP4 = -effCPD; }
    else if ( process==DDBarSSBar ){ cP3 = effCPS; cP4 = -effCPS; }
    else if ( process==SSBarUUBar ){ cP3 = effCPU; cP4 = -effCPU; }
    else if ( process==SSBarDDBar ){ cP3 = effCPD; cP4 = -effCPD; }
    else if ( process==SSBarSSBar ){ cP3 = effCPS; cP4 = -effCPS; }
    else if ( process==UBarUBarUBarUBar ){ cP3 = -effCPU; cP4 = -effCPU; }
    else if ( process==DBarDBarDBarDBar ){ cP3 = -effCPD; cP4 = -effCPD; }
    else if ( process==SBarSBarSBarSBar ){ cP3 = -effCPS; cP4 = -effCPS; }
    else if ( process==UBarDBarUBarDBar ){ cP3 = -effCPU; cP4 = -effCPD; }
    else if ( process==UBarSBarUBarSBar ){ cP3 = -effCPU; cP4 = -effCPS; }
}


double xij(double x, double eta, double mi, double mj)
{
    double xijAux = x - pow(mi + eta*mj,2);

    return xijAux;
}


////////////////////////////////////////////////////////////////////////////////////////
//s channel


double sChannelZeroMomentum(double s)
{
    double k0 = sqrt(s);

    return k0;
}


double sChannelThreeMomentum()
{   
    double k = 0.0;
    
    return k;
}


////////////////////////////////////////////////////////////////////////////////////////
//t channel


double tChannelZeroMomentum(double s, double m1, double m2, double m3, double m4)
{
    double k0 = ( pow(m1,2) - pow(m2,2) - pow(m3,2) + pow(m4,2) )/( 2.0*sqrt(s) );

    return k0;
}


double tChannelThreeMomentum(double s, double t, double m1, double m2, double m3, double m4)
{   
    double k = sqrt( fabs(pow(pow(m1,2) - pow(m2,2) - pow(m3,2) + pow(m4,2), 2)/( 4.0*s ) - t) );
    
    return k;
}


////////////////////////////////////////////////////////////////////////////////////////
//u channel


double uChannel(double s, double t, double m1, double m2, double m3, double m4)
{
    double uChannelAux = pow(m1,2) + pow(m2,2) + pow(m3,2) + pow(m4,2) - s - t;

    return uChannelAux;
}


double uChannelZeroMomentum(double s, double m1, double m2, double m3, double m4)
{
    double k0 = ( pow(m1,2) - pow(m2,2) + pow(m3,2) - pow(m4,2) )/( 2.0*sqrt(s) );

    return k0;
}


double uChannelThreeMomentum(double s, double t, double m1, double m2, double m3, double m4)
{   
    double k = sqrt( fabs(pow(pow(m1,2) - pow(m2,2) + pow(m3,2) - pow(m4,2), 2)/( 4.0*s ) - uChannel(s, t, m1, m2, m3, m4) ) );
    
    return k;
}


////////////////////////////////////////////////////////////////////////////////////////
//scattering angle


//known as triangle function lambda(x,y,z)=x^2+y^2+z^2-2xy-2xz-2yz
double lambdaTriangle(double x, double y, double z)
{
    double lambda = pow(x,2) + pow(y,2) + pow(z,2) - 2*x*y - 2*x*z - 2*y*z;

    return lambda;
}


double cosScatteringAngle(double s, double t, double m1, double m2, double m3, double m4)
{   
    double u = uChannel(s, t, m1, m2, m3, m4);

    double lambda12 = lambdaTriangle(s, pow(m1,2), pow(m2,2));

    double lambda34 = lambdaTriangle(s, pow(m3,2), pow(m4,2));

    double cosTheta = ( s*( t - u ) + ( pow(m1,2) - pow(m2,2) )*( pow(m3,2) - pow(m4,2) ) )/( sqrt(lambda12)*sqrt(lambda34) );

    return cosTheta;
}


double sinScatteringAngle(double s, double t, double m1, double m2, double m3, double m4)
{       
    double cosTheta = cosScatteringAngle(s, t, m1, m2, m3, m4);
    double sqrtArg = 1.0 - pow(cosTheta, 2);

    if ( sqrtArg<0 && fabs(sqrtArg)>1E-15 )
    { 
        cout << "Argument of the square root in sinScatteringAngle is negative and its modulus is bigger then 1E-15!\n"; 
    }

    double sinTheta = sqrt( fabs( sqrtArg ) );

    return sinTheta;
}


////////////////////////////////////////////////////////////////////////////////////////
//Calculate the sum of matrix elements involved in the quark-quark scatteing: 1/(4Nc^2) SUM |Mu - Mt|^2
double dsigmadtNJLQuarkQuarkScattering(double Nc, double s, double t, double m1, double m2, double m3, double m4, gsl_complex DPu, gsl_complex DSu, gsl_complex DPt, gsl_complex DSt)
{
    //calculate u channel from s and t channels
    double u = uChannel(s, t, m1, m2, m3, m4);

    ////////////////////////////////////////
    //uu channel term
    // 1/(4Nc^2) SUM |Muu|^2 = |DSu|^2 (u+14) (u+23) + |DPu|^2 (u−14) (u−23)

    gsl_complex auxMuu_1;
    double Xuu_1 = xij(u, +1, m1, m4)*xij(u, +1, m2, m3);
    auxMuu_1 = gsl_complex_mul(DSu, gsl_complex_conjugate(DSu));
    auxMuu_1 = gsl_complex_mul_real(auxMuu_1, +Xuu_1);

    gsl_complex auxMuu_2;
    double Xuu_2 = xij(u, -1, m1, m4)*xij(u, -1, m2, m3);
    auxMuu_2 = gsl_complex_mul(DPu, gsl_complex_conjugate(DPu));
    auxMuu_2 = gsl_complex_mul_real(auxMuu_2, +Xuu_2);

    gsl_complex auxMuu;
    auxMuu = gsl_complex_add(auxMuu_1, auxMuu_2);

    ////////////////////////////////////////
    //ut (interference) channel term
    // 1/(4Nc^2) SUM |Mut|^2 = 1/(4Nc) {  DSt DSu* [ (t+13) (t+24) - (s+12) (s+34) + (u+14) (u+23) ] 
    //                                  - DSt DPu* [ (t+13) (t+24) - (s-12) (s-34) + (u-14) (u-23) ] 
    //                                  - DPt DSu* [ (t-13) (t-24) - (s-12) (s-34) + (u+14) (u+23) ] 
    //                                  + DPt DPu* [ (t-13) (t-24) - (s+12) (s+34) + (u-14) (u-23) ]  }

    gsl_complex auxMut_1;
    double Xut_1 = xij(t, +1, m1, m3)*xij(t, +1, m2, m4) - xij(s, +1, m1, m2)*xij(s, +1, m3, m4) + xij(u, +1, m1, m4)*xij(u, +1, m2, m3);
    auxMut_1 = gsl_complex_mul(DSt, gsl_complex_conjugate(DSu));
    auxMut_1 = gsl_complex_mul_real(auxMut_1, +Xut_1/(4.0*Nc) );

    gsl_complex auxMut_2;
    double Xut_2 = xij(t, +1, m1, m3)*xij(t, +1, m2, m4) - xij(s, -1, m1, m2)*xij(s, -1, m3, m4) + xij(u, -1, m1, m4)*xij(u, -1, m2, m3);
    auxMut_2 = gsl_complex_mul(DSt, gsl_complex_conjugate(DPu));
    auxMut_2 = gsl_complex_mul_real(auxMut_2, -Xut_2/(4.0*Nc) );

    gsl_complex auxMut_3;
    double Xut_3 = xij(t, -1, m1, m3)*xij(t, -1, m2, m4) - xij(s, -1, m1, m2)*xij(s, -1, m3, m4) + xij(u, +1, m1, m4)*xij(u, +1, m2, m3);
    auxMut_3 = gsl_complex_mul(DPt, gsl_complex_conjugate(DSu));
    auxMut_3 = gsl_complex_mul_real(auxMut_3, -Xut_3/(4.0*Nc) );

    gsl_complex auxMut_4;
    double Xut_4 = xij(t, -1, m1, m3)*xij(t, -1, m2, m4) - xij(s, +1, m1, m2)*xij(s, +1, m3, m4) + xij(u, -1, m1, m4)*xij(u, -1, m2, m3);
    auxMut_4 = gsl_complex_mul(DPt, gsl_complex_conjugate(DPu));
    auxMut_4 = gsl_complex_mul_real(auxMut_4, +Xut_4/(4.0*Nc) );

    gsl_complex auxMut;
    auxMut = gsl_complex_add(auxMut_1, auxMut_2);
    auxMut = gsl_complex_add(auxMut, auxMut_3);
    auxMut = gsl_complex_add(auxMut, auxMut_4);

    ////////////////////////////////////////
    //tt channel term
    // 1/(4Nc^2) SUM |Mtt|^2 = |DSt|^2 t+13 t+24 + |DPt|^2 t−13 t−24

    gsl_complex auxMtt_1;
    double Xtt_1 = xij(t, +1, m1, m3)*xij(t, +1, m2, m4);
    auxMtt_1 = gsl_complex_mul(DSt, gsl_complex_conjugate(DSt));
    auxMtt_1 = gsl_complex_mul_real(auxMtt_1, Xtt_1);

    gsl_complex auxMtt_2;
    double Xtt_2 = xij(t, -1, m1, m3)*xij(t, -1, m2, m4);
    auxMtt_2 = gsl_complex_mul(DPt, gsl_complex_conjugate(DPt));
    auxMtt_2 = gsl_complex_mul_real(auxMtt_2, Xtt_2);

    gsl_complex auxMtt = gsl_complex_add(auxMtt_1, auxMtt_2);
    
    ////////////////////////////////////////
    // 1/(4Nc^2) SUM |Mu - Mt|^2 = |Muu|^2 + |Mtt|^2 - 2Re[Mut]
    //(z1-z2)(z1-z2)* = (z1-z2)(z1*-z2*) = z1 z1* - z1 z2* - z2z1* + z2z2* = z1 z1* + z2z2* - 2Re[z1 z2*]
    double auxM = GSL_REAL(auxMuu) + GSL_REAL(auxMtt) - 2.0*GSL_REAL(auxMut);

    // dsigma/dt = 1/(16Pi^2 (s+12) (s-12)) 1/(4Nc^2) SUM |Mu - Mt|^2
    double dsigmadt = auxM/( 16.0*M_PI*xij(s, +1, m1, m2)*xij(s, -1, m1, m2) );

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//Calculate the sum of matrix elements involved in the quark-antiquark scatteing: 1/(4Nc^2) SUM |Ms - Mt|^2
double dsigmadtNJLQuarkAntiquarkScattering(double Nc, double s, double t, double m1, double m2, double m3, double m4, gsl_complex DPs, gsl_complex DSs, gsl_complex DPt, gsl_complex DSt)
{
    //calculate u channel from s and t channels
    double u = uChannel(s, t, m1, m2, m3, m4);

    ////////////////////////////////////////
    //ss channel term
    // 1/(4Nc^2) SUM |Mss|^2 = |DSs|^2 (s+12) (s+34) + |DPs|^2 (s−12) (s−34)

    gsl_complex auxMss_1;
    double Xss_1 = xij(s, +1, m1, m2)*xij(s, +1, m3, m4);
    auxMss_1 = gsl_complex_mul(DSs, gsl_complex_conjugate(DSs));
    auxMss_1 = gsl_complex_mul_real(auxMss_1, +Xss_1);

    gsl_complex auxMss_2;
    double Xss_2 = xij(s, -1, m1, m2)*xij(s, -1, m3, m4);
    auxMss_2 = gsl_complex_mul(DPs, gsl_complex_conjugate(DPs));
    auxMss_2 = gsl_complex_mul_real(auxMss_2, +Xss_2);

    gsl_complex auxMss;
    auxMss = gsl_complex_add(auxMss_1, auxMss_2);

    ////////////////////////////////////////
    //ut (interference) channel term
    // 1/(4Nc^2) SUM |Mst|^2 = 1/(4Nc) {  DSs DSt* [ (s+12) (s+34) − (u+14) (u+23) + (t+13) (t+24) ]
    //                                  − DSs DPt* [ (s+12) (s+34) − (u−14) (u−23) + (t−13) (t−24) ]
    //                                  − DPs DSt* [ (s−12) (s−34) − (u−14) (u−23) + (t+13) (t+24) ] 
    //                                  + DPs DPt* [ (s−12) (s−34) − (u+14) (u+23) + (t−13) (t−24) ]  }

    gsl_complex auxMst_1;
    double Xst_1 = xij(s, +1, m1, m2)*xij(s, +1, m3, m4) - xij(u, +1, m1, m4)*xij(u, +1, m2, m3) + xij(t, +1, m1, m3)*xij(t, +1, m2, m4);
    auxMst_1 = gsl_complex_mul(DSs, gsl_complex_conjugate(DSt));
    auxMst_1 = gsl_complex_mul_real(auxMst_1, +Xst_1/(4.0*Nc) );

    gsl_complex auxMst_2;
    double Xst_2 = xij(s, +1, m1, m2)*xij(s, +1, m3, m4) - xij(u, -1, m1, m4)*xij(u, -1, m2, m3) + xij(t, -1, m1, m3)*xij(t, -1, m2, m4);
    auxMst_2 = gsl_complex_mul(DSs, gsl_complex_conjugate(DPt));
    auxMst_2 = gsl_complex_mul_real(auxMst_2, -Xst_2/(4.0*Nc) );

    gsl_complex auxMst_3;
    double Xst_3 = xij(s, -1, m1, m2)*xij(s, -1, m3, m4) - xij(u, -1, m1, m4)*xij(u, -1, m2, m3) + xij(t, +1, m1, m3)*xij(t, +1, m2, m4);
    auxMst_3 = gsl_complex_mul(DPs, gsl_complex_conjugate(DSt));
    auxMst_3 = gsl_complex_mul_real(auxMst_3, -Xst_3/(4.0*Nc) );

    gsl_complex auxMst_4;
    double Xst_4 = xij(s, -1, m1, m2)*xij(s, -1, m3, m4) - xij(u, +1, m1, m4)*xij(u, +1, m2, m3) + xij(t, -1, m1, m3)*xij(t, -1, m2, m4);
    auxMst_4 = gsl_complex_mul(DPs, gsl_complex_conjugate(DPt));
    auxMst_4 = gsl_complex_mul_real(auxMst_4, +Xst_4/(4.0*Nc) );

    gsl_complex auxMst;
    auxMst = gsl_complex_add(auxMst_1, auxMst_2);
    auxMst = gsl_complex_add(auxMst, auxMst_3);
    auxMst = gsl_complex_add(auxMst, auxMst_4);

    ////////////////////////////////////////
    //tt channel term
    // 1/(4Nc^2) SUM |Mtt|^2 = |DSt|^2 (t+13) (t+24) + |DPt|^2 (t−13) (t−24)

    gsl_complex auxMtt_1;
    double Xtt_1 = xij(t, +1, m1, m3)*xij(t, +1, m2, m4);
    auxMtt_1 = gsl_complex_mul(DSt, gsl_complex_conjugate(DSt));
    auxMtt_1 = gsl_complex_mul_real(auxMtt_1, Xtt_1);

    gsl_complex auxMtt_2;
    double Xtt_2 = xij(t, -1, m1, m3)*xij(t, -1, m2, m4);
    auxMtt_2 = gsl_complex_mul(DPt, gsl_complex_conjugate(DPt));
    auxMtt_2 = gsl_complex_mul_real(auxMtt_2, Xtt_2);

    gsl_complex auxMtt = gsl_complex_add(auxMtt_1, auxMtt_2);

    ////////////////////////////////////////
    // 1/(4Nc^2) SUM |Ms - Mt|^2 = |Mss|^2 + |Mtt|^2 - 2Re[Mst]
    //(z1-z2)(z1-z2)* = (z1-z2)(z1*-z2*) = z1 z1* - z1 z2* - z2z1* + z2z2* = z1 z1* + z2z2* - 2Re[z1 z2*]
    double auxM = GSL_REAL(auxMss) + GSL_REAL(auxMtt) - 2.0*GSL_REAL(auxMst);

    // dsigma/dt = 1/(16Pi^2 (s+12) (s-12)) 1/(4Nc^2) SUM |Ms - Mt|^2
    double dsigmadt = auxM/( 16.0*M_PI*xij(s, +1, m1, m2)*xij(s, -1, m1, m2) );

    return dsigmadt;
}



////////////////////////////////////////////////////////////////////////////////////////
//UD->UD differential cross section
double dsigmadtUDUD(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassU;
    double m2 = effMassD;
    double m3 = effMassU;
    double m4 = effMassD;

    //calculate propagators that take part in this process
    double uK0 = uChannelZeroMomentum(s, m1, m2, m3, m4);
    double uKVec = uChannelThreeMomentum(s, t, m1, m2, m3, m4);
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);


    // DPu = 2*pionPlusPropagator
    gsl_complex DPu;
    DPu = pionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DPu = gsl_complex_mul_real(DPu, 2);


    // DSu = 2 Dsigmapionplus
    gsl_complex DSu;
    DSu = sigmaPionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DSu = gsl_complex_mul_real(DSu, 2);


    // DPt = (2./3.)*DnP00 + (2sqrt(2)/3)*DnP08 - DnP33 + (1./3.)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    //gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    gsl_complex DnPt33 = DnPt.getValue(1, 1);
    //gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, +2.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt33, -1.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +1.0/3.0));
    

    // DSt = (2./3.)*DnS00 + (2sqrt(2)/3)*DnS08 - DnS33 + (1./3.)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    //gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    gsl_complex DnSt33 = DnSt.getValue(1, 1);
    //gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, +2.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt33, -1.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, 1.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkQuarkScattering(Nc, s, t, m1, m2, m3, m4, DPu, DSu, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//DU->DU differential cross section
double dsigmadtDUDU(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassD;
    double m2 = effMassU;
    double m3 = effMassD;
    double m4 = effMassU;

    //calculate propagators that take part in this process
    double uK0 = uChannelZeroMomentum(s, m1, m2, m3, m4);
    double uKVec = uChannelThreeMomentum(s, t, m1, m2, m3, m4);
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPu = 2*pionMinusPropagator
    gsl_complex DPu;
    DPu = pionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DPu = gsl_complex_mul_real(DPu, 2);


    // DSu = 2*sigmaPionMinusPropagator
    gsl_complex DSu;
    DSu = sigmaPionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DSu = gsl_complex_mul_real(DSu, 2);


    // DPt = (2./3.)*DnP00 + (2sqrt(2)/3)*DnP08 - DnP33 + (1./3.)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    //gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    gsl_complex DnPt33 = DnPt.getValue(1, 1);
    //gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, +2.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt33, -1.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +1.0/3.0));
    

    // DSt = (2./3.)*DnS00 + (2sqrt(2)/3)*DnS08 - DnS33 + (1./3.)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    //gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    gsl_complex DnSt33 = DnSt.getValue(1, 1);
    //gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, +2.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt33, -1.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, 1.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkQuarkScattering(Nc, s, t, m1, m2, m3, m4, DPu, DSu, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//US->US differential cross section
double dsigmadtUSUS(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassU;
    double m2 = effMassS;
    double m3 = effMassU;
    double m4 = effMassS;

    //calculate propagators that take part in this process
    double uK0 = uChannelZeroMomentum(s, m1, m2, m3, m4);
    double uKVec = uChannelThreeMomentum(s, t, m1, m2, m3, m4);
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPu = 2*kaonPlusPropagator
    gsl_complex DPu;
    DPu = kaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DPu = gsl_complex_mul_real(DPu, 2);


    // DSu = 2*sigmaKaonPlusPropagator
    gsl_complex DSu;
    DSu = sigmaKaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DSu = gsl_complex_mul_real(DSu, 2);


    // DPt = (2./3.)*DnP00 + sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 - (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, +sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, -2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, -2.0/3.0));
    

    // DSt = (2./3.)*DnS00 + sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 - (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, +sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, -2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, -2.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkQuarkScattering(Nc, s, t, m1, m2, m3, m4, DPu, DSu, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//SU->SU differential cross section
double dsigmadtSUSU(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassS;
    double m2 = effMassU;
    double m3 = effMassS;
    double m4 = effMassU;

    //calculate propagators that take part in this process
    double uK0 = uChannelZeroMomentum(s, m1, m2, m3, m4);
    double uKVec = uChannelThreeMomentum(s, t, m1, m2, m3, m4);
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPu = 2*kaonMinusPropagator
    gsl_complex DPu;
    DPu = kaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DPu = gsl_complex_mul_real(DPu, 2);


    // DSu = 2*sigmaKaonMinusPropagator
    gsl_complex DSu;
    DSu = sigmaKaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DSu = gsl_complex_mul_real(DSu, 2);


    // DPt = (2./3.)*DnP00 + sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 - (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, +sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, -2./sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, -2.0/3.0));
    

    // DSt = (2./3.)*DnS00 + sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 - (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, +sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, -2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, -2.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkQuarkScattering(Nc, s, t, m1, m2, m3, m4, DPu, DSu, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//DS->DS differential cross section
double dsigmadtDSDS(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassD;
    double m2 = effMassS;
    double m3 = effMassD;
    double m4 = effMassS;

    //calculate propagators that take part in this process
    double uK0 = uChannelZeroMomentum(s, m1, m2, m3, m4);
    double uKVec = uChannelThreeMomentum(s, t, m1, m2, m3, m4);
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPu = 2*neutralKaonPropagator
    gsl_complex DPu;
    DPu = neutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DPu = gsl_complex_mul_real(DPu, 2);


    // DSu = 2*neutralSigmaKaonPropagator
    gsl_complex DSu;
    DSu = neutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DSu = gsl_complex_mul_real(DSu, 2);


    // DPt = (2./3.)*DnP00 - sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 + (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, -sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, +2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, -2.0/3.0));
    

    // DSt = (2./3.)*DnS00 - sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 + (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, -sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, +2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, -2.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkQuarkScattering(Nc, s, t, m1, m2, m3, m4, DPu, DSu, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//SD->SD differential cross section
double dsigmadtSDSD(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassS;
    double m2 = effMassD;
    double m3 = effMassS;
    double m4 = effMassD;

    //calculate propagators that take part in this process
    double uK0 = uChannelZeroMomentum(s, m1, m2, m3, m4);
    double uKVec = uChannelThreeMomentum(s, t, m1, m2, m3, m4);
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPu = 2*antiNeutralKaonPropagator
    gsl_complex DPu;
    DPu = antiNeutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DPu = gsl_complex_mul_real(DPu, 2);


    // DSu = 2*antiNeutralSigmaKaonPropagator
    gsl_complex DSu;
    DSu = antiNeutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    DSu = gsl_complex_mul_real(DSu, 2);


    // DPt = (2./3.)*DnP00 - sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 + (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, -sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, +2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, -2.0/3.0));
    

    // DSt = (2./3.)*DnS00 - sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 + (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, -sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, +2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, -2.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkQuarkScattering(Nc, s, t, m1, m2, m3, m4, DPu, DSu, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//UU->UU differential cross section
double dsigmadtUUUU(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassU;
    double m2 = effMassU;
    double m3 = effMassU;
    double m4 = effMassU;

    //calculate propagators that take part in this process
    double uK0 = uChannelZeroMomentum(s, m1, m2, m3, m4);
    double uKVec = uChannelThreeMomentum(s, t, m1, m2, m3, m4);
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPu = (2./3.)*DnP00 + 2*sqrt(2/3)*DnP03 + (2*sqrt(2)/3)*DnP08 + DnP33 + (2/sqrt(3))*DnP38 + (1/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPu = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision); 
    gsl_complex DnPu00 = DnPu.getValue(0, 0);
    gsl_complex DnPu03 = DnPu.getValue(0, 1);
    gsl_complex DnPu08 = DnPu.getValue(0, 2);
    gsl_complex DnPu33 = DnPu.getValue(1, 1);
    gsl_complex DnPu38 = DnPu.getValue(1, 2);
    gsl_complex DnPu88 = DnPu.getValue(2, 2);

    gsl_complex DPu;
    DPu = gsl_complex_mul_real(DnPu00, +2.0/3.0);
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu03, +2.0*sqrt(2.0/3.0)));
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu08, +2.0*sqrt(2.0)/3.0));
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu33, +1.0));
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu38, +2.0/sqrt(3.0)));
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu88, +1.0/3.0));


    // DSu = (2./3.)*DnS00 + 2*sqrt(2/3)*DnS03 + (2*sqrt(2)/3)*DnS08 + DnS33 + (2/sqrt(3))*DnS38 + (1/3)*DnS88 ,  where DnS is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSu = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    gsl_complex DnSu00 = DnSu.getValue(0, 0);
    gsl_complex DnSu03 = DnSu.getValue(0, 1);
    gsl_complex DnSu08 = DnSu.getValue(0, 2);
    gsl_complex DnSu33 = DnSu.getValue(1, 1);
    gsl_complex DnSu38 = DnSu.getValue(1, 2);
    gsl_complex DnSu88 = DnSu.getValue(2, 2);

    gsl_complex DSu;
    DSu = gsl_complex_mul_real(DnSu00, +2.0/3.0);
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu03, +2.0*sqrt(2.0/3.0)));
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu08, +2.0*sqrt(2.0)/3.0));
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu33, +1.0));
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu38, +2.0/sqrt(3.0)));
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu88, +1.0/3.0));


    // DPt = (2./3.)*DnP00 + 2*sqrt(2/3)*DnP03 + (2*sqrt(2)/3)*DnP08 + DnP33 + (2/sqrt(3))*DnP38 + (1/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, +2.0*sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, +2.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt33, +1.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, +2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +1.0/3.0));

    
    // DSt = (2./3.)*DnS00 + 2*sqrt(2/3)*DnS03 + (2*sqrt(2)/3)*DnS08 + DnS33 + (2/sqrt(3))*DnS38 + (1/3)*DnS88 ,  where DnP is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, +2.0*sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, +2.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt33, +1.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, +2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, +1.0/3.0));


    ////calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkQuarkScattering(Nc, s, t, m1, m2, m3, m4, DPu, DSu, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//DD->DD differential cross section
double dsigmadtDDDD(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassD;
    double m2 = effMassD;
    double m3 = effMassD;
    double m4 = effMassD;

    //calculate propagators that take part in this process
    double uK0 = uChannelZeroMomentum(s, m1, m2, m3, m4);
    double uKVec = uChannelThreeMomentum(s, t, m1, m2, m3, m4);
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPu = (2./3.)*DnP00 - 2*sqrt(2/3)*DnP03 + (2*sqrt(2)/3)*DnP08 + DnP33 - (2/sqrt(3))*DnP38 + (1/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPu = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision); 
    gsl_complex DnPu00 = DnPu.getValue(0, 0);
    gsl_complex DnPu03 = DnPu.getValue(0, 1);
    gsl_complex DnPu08 = DnPu.getValue(0, 2);
    gsl_complex DnPu33 = DnPu.getValue(1, 1);
    gsl_complex DnPu38 = DnPu.getValue(1, 2);
    gsl_complex DnPu88 = DnPu.getValue(2, 2);

    gsl_complex DPu;
    DPu = gsl_complex_mul_real(DnPu00, +2.0/3.0);
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu03, -2.0*sqrt(2.0/3.0)));
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu08, +2.0*sqrt(2.0)/3.0));
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu33, +1.0));
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu38, -2.0/sqrt(3.0)));
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu88, +1.0/3.0));


    // DSu = (2./3.)*DnS00 - 2*sqrt(2/3)*DnS03 + (2*sqrt(2)/3)*DnS08 + DnS33 - (2/sqrt(3))*DnS38 + (1/3)*DnS88 ,  where DnS is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSu = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    gsl_complex DnSu00 = DnSu.getValue(0, 0);
    gsl_complex DnSu03 = DnSu.getValue(0, 1);
    gsl_complex DnSu08 = DnSu.getValue(0, 2);
    gsl_complex DnSu33 = DnSu.getValue(1, 1);
    gsl_complex DnSu38 = DnSu.getValue(1, 2);
    gsl_complex DnSu88 = DnSu.getValue(2, 2);

    gsl_complex DSu;
    DSu = gsl_complex_mul_real(DnSu00, +2.0/3.0);
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu03, -2.0*sqrt(2.0/3.0)));
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu08, +2.0*sqrt(2.0)/3.0));
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu33, +1.0));
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu38, -2.0/sqrt(3.0)));
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu88, +1.0/3.0));


    // DPt = (2./3.)*DnP00 - 2*sqrt(2/3)*DnP03 + (2*sqrt(2)/3)*DnP08 + DnP33 - (2/sqrt(3))*DnP38 + (1/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, -2.0*sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, +2.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt33, +1.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, -2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +1.0/3.0));

    
    // DSt = (2./3.)*DnS00 - 2*sqrt(2/3)*DnS03 + (2*sqrt(2)/3)*DnS08 + DnS33 - (2/sqrt(3))*DnS38 + (1/3)*DnS88 ,  where DnP is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, -2.0*sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, +2.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt33, +1.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, -2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, +1.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkQuarkScattering(Nc, s, t, m1, m2, m3, m4, DPu, DSu, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//SS->SS differential cross section
double dsigmadtSSSS(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassS;
    double m2 = effMassS;
    double m3 = effMassS;
    double m4 = effMassS;

    //calculate propagators that take part in this process
    double uK0 = uChannelZeroMomentum(s, m1, m2, m3, m4);
    double uKVec = uChannelThreeMomentum(s, t, m1, m2, m3, m4);
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPu = (2./3.)*DnP00 - (4*sqrt(2)/3)*DnP08 + (4/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPu = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision); 
    gsl_complex DnPu00 = DnPu.getValue(0, 0);
    //gsl_complex DnPu03 = DnPu.getValue(0, 1);
    gsl_complex DnPu08 = DnPu.getValue(0, 2);
    //gsl_complex DnPu33 = DnPu.getValue(1, 1);
    //gsl_complex DnPu38 = DnPu.getValue(1, 2);
    gsl_complex DnPu88 = DnPu.getValue(2, 2);

    gsl_complex DPu;
    DPu = gsl_complex_mul_real(DnPu00, +2.0/3.0);
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu08, -4.0*sqrt(2.0)/3.0));
    DPu = gsl_complex_add(DPu, gsl_complex_mul_real(DnPu88, +4.0/3.0));


    // DSu = (2./3.)*DnS00 - (4*sqrt(2)/3)*DnS08 + (4/3)*DnS88 ,  where DnS is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSu = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, uK0, uKVec, integralPrecision);
    gsl_complex DnSu00 = DnSu.getValue(0, 0);
    //gsl_complex DnSu03 = DnSu.getValue(0, 1);
    gsl_complex DnSu08 = DnSu.getValue(0, 2);
    //gsl_complex DnSu33 = DnSu.getValue(1, 1);
    //gsl_complex DnSu38 = DnSu.getValue(1, 2);
    gsl_complex DnSu88 = DnSu.getValue(2, 2);

    gsl_complex DSu;
    DSu = gsl_complex_mul_real(DnSu00, +2.0/3.0);
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu08, -4.0*sqrt(2.0)/3.0));
    DSu = gsl_complex_add(DSu, gsl_complex_mul_real(DnSu88, +4.0/3.0));


    // DPt = (2./3.)*DnP00 - (4*sqrt(2)/3)*DnP08 + (4/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    //gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    //gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -4.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +4.0/3.0));

    
    // DSt = (2./3.)*DnS00 - (4*sqrt(2)/3)*DnS08 + (4/3)*DnS88 ,  where DnP is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    //gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    //gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -4.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, +4.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkQuarkScattering(Nc, s, t, m1, m2, m3, m4, DPu, DSu, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//UDBar->UDBar differential cross section
double dsigmadtUDBarUDBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassU;
    double m2 = effMassD;
    double m3 = effMassU;
    double m4 = effMassD;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = 2*pionPlusPropagator
    gsl_complex DPs;
    DPs = pionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DPs = gsl_complex_mul_real(DPs, 2);


    // DSs = 2*sigmaPionPlusPropagator
    gsl_complex DSs;
    DSs = sigmaPionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DSs = gsl_complex_mul_real(DSs, 2);


    // DPt = (2./3.)*DnP00 + (2sqrt(2)/3)*DnP08 - DnP33 + (1./3.)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    //gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    gsl_complex DnPt33 = DnPt.getValue(1, 1);
    //gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, +2.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt33, -1.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +1.0/3.0));
    

    // DSt = (2./3.)*DnS00 + (2sqrt(2)/3)*DnS08 - DnS33 + (1./3.)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    //gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    gsl_complex DnSt33 = DnSt.getValue(1, 1);
    //gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, +2.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt33, -1.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, 1.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//DUBar->DUBar differential cross section
double dsigmadtDUBarDUBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassD;
    double m2 = effMassU;
    double m3 = effMassD;
    double m4 = effMassU;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = 2*pionMinusPropagator
    gsl_complex DPs;
    DPs = pionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DPs = gsl_complex_mul_real(DPs, 2);


    // DSs = 2*sigmaPionMinusPropagator
    gsl_complex DSs;
    DSs = sigmaPionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DSs = gsl_complex_mul_real(DSs, 2);


    // DPt = (2./3.)*DnP00 + (2sqrt(2)/3)*DnP08 - DnP33 + (1./3.)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    //gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    gsl_complex DnPt33 = DnPt.getValue(1, 1);
    //gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, +2.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt33, -1.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +1.0/3.0));
    

    // DSt = (2./3.)*DnS00 + (2sqrt(2)/3)*DnS08 - DnS33 + (1./3.)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    //gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    gsl_complex DnSt33 = DnSt.getValue(1, 1);
    //gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, +2.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt33, -1.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, 1.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//USBar->USBar differential cross section
double dsigmadtUSBarUSBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassU;
    double m2 = effMassS;
    double m3 = effMassU;
    double m4 = effMassS;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = 2*kaonPlusPropagator
    gsl_complex DPs;
    DPs = kaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DPs = gsl_complex_mul_real(DPs, 2);


    // DSs = 2*sigmaKaonPlusPropagator
    gsl_complex DSs;
    DSs = sigmaKaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DSs = gsl_complex_mul_real(DSs, 2);


    // DPt = (2./3.)*DnP00 + sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 - (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, +sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, -2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, -2.0/3.0));
    

    // DSt = (2./3.)*DnS00 + sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 - (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, +sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, -2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, -2.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//SUBar->SUBar differential cross section
double dsigmadtSUBarSUBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassS;
    double m2 = effMassU;
    double m3 = effMassS;
    double m4 = effMassU;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = 2*kaonMinusPropagator
    gsl_complex DPs;
    DPs = kaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DPs = gsl_complex_mul_real(DPs, 2);


    // DSs = 2*sigmaKaonMinusPropagator
    gsl_complex DSs;
    DSs = sigmaKaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DSs = gsl_complex_mul_real(DSs, 2);


    // DPt = (2./3.)*DnP00 + sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 - (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, +sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, -2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, -2.0/3.0));
    

    // DSt = (2./3.)*DnS00 + sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 - (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, +sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, -2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, -2.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//DSBar->DSBar differential cross section
double dsigmadtDSBarDSBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                    double effChemPotU, double effChemPotD, double effChemPotS, 
                    double effMassU, double effMassD, double effMassS, 
                    double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassD;
    double m2 = effMassS;
    double m3 = effMassD;
    double m4 = effMassS;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = 2*neutralKaonPropagator
    gsl_complex DPs;
    DPs = neutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DPs = gsl_complex_mul_real(DPs, 2);


    // DSs = 2*neutralSigmaKaonPropagator
    gsl_complex DSs;
    DSs = neutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DSs = gsl_complex_mul_real(DSs, 2);


    // DPt = (2./3.)*DnP00 - sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 + (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, -sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, +2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, -2.0/3.0));
    

    // DSt = (2./3.)*DnS00 - sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 + (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, -sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, +2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, -2.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//SDBar->SDBar differential cross section
double dsigmadtSDBarSDBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassS;
    double m2 = effMassD;
    double m3 = effMassS;
    double m4 = effMassD;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = 2*antiNeutralKaonPropagator
    gsl_complex DPs;
    DPs = antiNeutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DPs = gsl_complex_mul_real(DPs, 2);


    // DSs = 2*antiNeutralSigmaKaonPropagator
    gsl_complex DSs;
    DSs = antiNeutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    DSs = gsl_complex_mul_real(DSs, 2);


    // DPt = (2./3.)*DnP00 - sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 + (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, -sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, +2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, -2.0/3.0));
    

    // DSt = (2./3.)*DnS00 - sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 + (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    // DSt = (2./3.)*DnS00 - sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 + (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, -sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, +2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, -2.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//UUBar->UUBar differential cross section
double dsigmadtUUBarUUBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassU;
    double m2 = effMassU;
    double m3 = effMassU;
    double m4 = effMassU;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = (2./3.)*DnP00 + 2*sqrt(2/3)*DnP03 + (2*sqrt(2)/3)*DnP08 + DnP33 + (2/sqrt(3))*DnP38 + (1/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPs = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision); 
    gsl_complex DnPs00 = DnPs.getValue(0, 0);
    gsl_complex DnPs03 = DnPs.getValue(0, 1);
    gsl_complex DnPs08 = DnPs.getValue(0, 2);
    gsl_complex DnPs33 = DnPs.getValue(1, 1);
    gsl_complex DnPs38 = DnPs.getValue(1, 2);
    gsl_complex DnPs88 = DnPs.getValue(2, 2);

    gsl_complex DPs;
    DPs = gsl_complex_mul_real(DnPs00, +2.0/3.0);
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs03, +2.0*sqrt(2.0/3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs08, +2.0*sqrt(2.0)/3.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs33, +1.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs38, +2.0/sqrt(3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs88, +1.0/3.0));


    // DSs = (2./3.)*DnS00 + 2*sqrt(2/3)*DnS03 + (2*sqrt(2)/3)*DnS08 + DnS33 + (2/sqrt(3))*DnS38 + (1/3)*DnS88 ,  where DnS is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSs = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    gsl_complex DnSs00 = DnSs.getValue(0, 0);
    gsl_complex DnSs03 = DnSs.getValue(0, 1);
    gsl_complex DnSs08 = DnSs.getValue(0, 2);
    gsl_complex DnSs33 = DnSs.getValue(1, 1);
    gsl_complex DnSs38 = DnSs.getValue(1, 2);
    gsl_complex DnSs88 = DnSs.getValue(2, 2);

    gsl_complex DSs;
    DSs = gsl_complex_mul_real(DnSs00, +2.0/3.0);
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs03, +2.0*sqrt(2.0/3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs08, +2.0*sqrt(2.0)/3.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs33, +1.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs38, +2.0/sqrt(3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs88, +1.0/3.0));


    // DPt = (2./3.)*DnP00 + 2*sqrt(2/3)*DnP03 + (2*sqrt(2)/3)*DnP08 + DnP33 + (2/sqrt(3))*DnP38 + (1/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, +2.0*sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, +2.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt33, +1.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, +2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +1.0/3.0));

    
    // DSt = (2./3.)*DnS00 + 2*sqrt(2/3)*DnS03 + (2*sqrt(2)/3)*DnS08 + DnS33 + (2/sqrt(3))*DnS38 + (1/3)*DnS88 ,  where DnP is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, +2.0*sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, +2.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt33, +1.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, +2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, +1.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//UUBar->DDBar differential cross section
double dsigmadtUUBarDDBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassU;
    double m2 = effMassU;
    double m3 = effMassD;
    double m4 = effMassD;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = (2./3.)*DnP00 + (2sqrt(2)/3)*DnP08 - DnP33 + (1./3.)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPs = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision); 
    gsl_complex DnPs00 = DnPs.getValue(0, 0);
    //gsl_complex DnPs03 = DnPs.getValue(0, 1);
    gsl_complex DnPs08 = DnPs.getValue(0, 2);
    gsl_complex DnPs33 = DnPs.getValue(1, 1);
    //gsl_complex DnPs38 = DnPs.getValue(1, 2);
    gsl_complex DnPs88 = DnPs.getValue(2, 2);

    gsl_complex DPs;
    DPs = gsl_complex_mul_real(DnPs00, +2.0/3.0);
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs08, +2.0*sqrt(2.0)/3.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs33, -1.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs88, +1.0/3.0));
    

    // DSs = (2./3.)*DnS00 + (2sqrt(2)/3)*DnS08 - DnS33 + (1./3.)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSs = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    gsl_complex DnSs00 = DnSs.getValue(0, 0);
    //gsl_complex DnSs03 = DnSs.getValue(0, 1);
    gsl_complex DnSs08 = DnSs.getValue(0, 2);
    gsl_complex DnSs33 = DnSs.getValue(1, 1);
    //gsl_complex DnSs38 = DnSs.getValue(1, 2);
    gsl_complex DnSs88 = DnSs.getValue(2, 2);

    gsl_complex DSs;
    DSs = gsl_complex_mul_real(DnSs00, +2.0/3.0);
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs08, +2.0*sqrt(2.0)/3.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs33, -1.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs88, 1.0/3.0));


    // DPt = 2*pionPlusPropagator
    gsl_complex DPt;
    DPt = pionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DPt = gsl_complex_mul_real(DPt, 2);


    // DSt = 2*sigmaPionPlusPropagator
    gsl_complex DSt;
    DSt = sigmaPionPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DSt = gsl_complex_mul_real(DSt, 2);


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//UUBar->SSBar differential cross section
double dsigmadtUUBarSSBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassU;
    double m2 = effMassU;
    double m3 = effMassS;
    double m4 = effMassS;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = (2./3.)*DnP00 + sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 - (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPs = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision); 
    gsl_complex DnPs00 = DnPs.getValue(0, 0);
    gsl_complex DnPs03 = DnPs.getValue(0, 1);
    gsl_complex DnPs08 = DnPs.getValue(0, 2);
    //gsl_complex DnPs33 = DnPs.getValue(1, 1);
    gsl_complex DnPs38 = DnPs.getValue(1, 2);
    gsl_complex DnPs88 = DnPs.getValue(2, 2);

    gsl_complex DPs;
    DPs = gsl_complex_mul_real(DnPs00, +2.0/3.0);
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs03, +sqrt(2.0/3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs08, -sqrt(2.0)/3.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs38, -2.0/sqrt(3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs88, -2.0/3.0));
    

    // DSs = (2./3.)*DnS00 + sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 - (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSs = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    gsl_complex DnSs00 = DnSs.getValue(0, 0);
    gsl_complex DnSs03 = DnSs.getValue(0, 1);
    gsl_complex DnSs08 = DnSs.getValue(0, 2);
    //gsl_complex DnSs33 = DnSs.getValue(1, 1);
    gsl_complex DnSs38 = DnSs.getValue(1, 2);
    gsl_complex DnSs88 = DnSs.getValue(2, 2);

    gsl_complex DSs;
    DSs = gsl_complex_mul_real(DnSs00, +2.0/3.0);
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs03, +sqrt(2.0/3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs08, -sqrt(2.0)/3.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs38, -2.0/sqrt(3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs88, -2.0/3.0));


    // DPt = 2*kaonPlusPropagator
    gsl_complex DPt;
    DPt = kaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DPt = gsl_complex_mul_real(DPt, 2);


    // DSt = 2*sigmaKaonPlusPropagator
    gsl_complex DSt;
    DSt = sigmaKaonPlusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DSt = gsl_complex_mul_real(DSt, 2);


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//DDBar->UUBar differential cross section
double dsigmadtDDBarUUBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassD;
    double m2 = effMassD;
    double m3 = effMassU;
    double m4 = effMassU;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = (2./3.)*DnP00 + (2sqrt(2)/3)*DnP08 - DnP33 + (1./3.)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPs = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision); 
    gsl_complex DnPs00 = DnPs.getValue(0, 0);
    //gsl_complex DnPs03 = DnPs.getValue(0, 1);
    gsl_complex DnPs08 = DnPs.getValue(0, 2);
    gsl_complex DnPs33 = DnPs.getValue(1, 1);
    //gsl_complex DnPs38 = DnPs.getValue(1, 2);
    gsl_complex DnPs88 = DnPs.getValue(2, 2);

    gsl_complex DPs;
    DPs = gsl_complex_mul_real(DnPs00, +2.0/3.0);
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs08, +2.0*sqrt(2.0)/3.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs33, -1.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs88, +1.0/3.0));
    

    // DSs = (2./3.)*DnS00 + (2sqrt(2)/3)*DnS08 - DnS33 + (1./3.)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSs = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    gsl_complex DnSs00 = DnSs.getValue(0, 0);
    //gsl_complex DnSs03 = DnSs.getValue(0, 1);
    gsl_complex DnSs08 = DnSs.getValue(0, 2);
    gsl_complex DnSs33 = DnSs.getValue(1, 1);
    //gsl_complex DnSs38 = DnSs.getValue(1, 2);
    gsl_complex DnSs88 = DnSs.getValue(2, 2);

    gsl_complex DSs;
    DSs = gsl_complex_mul_real(DnSs00, +2.0/3.0);
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs08, +2.0*sqrt(2.0)/3.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs33, -1.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs88, 1.0/3.0));


    // DPt = 2*pionMinusPropagator
    gsl_complex DPt;
    DPt = pionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DPt = gsl_complex_mul_real(DPt, 2);


    // DSt = 2*sigmaPionMinusPropagator
    gsl_complex DSt;
    DSt = sigmaPionMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DSt = gsl_complex_mul_real(DSt, 2);


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//DDBar->DDBar differential cross section
double dsigmadtDDBarDDBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassD;
    double m2 = effMassD;
    double m3 = effMassD;
    double m4 = effMassD;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = (2./3.)*DnP00 - 2*sqrt(2/3)*DnP03 + (2*sqrt(2)/3)*DnP08 + DnP33 - (2/sqrt(3))*DnP38 + (1/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPs = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision); 
    gsl_complex DnPs00 = DnPs.getValue(0, 0);
    gsl_complex DnPs03 = DnPs.getValue(0, 1);
    gsl_complex DnPs08 = DnPs.getValue(0, 2);
    gsl_complex DnPs33 = DnPs.getValue(1, 1);
    gsl_complex DnPs38 = DnPs.getValue(1, 2);
    gsl_complex DnPs88 = DnPs.getValue(2, 2);

    gsl_complex DPs;
    DPs = gsl_complex_mul_real(DnPs00, +2.0/3.0);
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs03, -2.0*sqrt(2.0/3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs08, +2.0*sqrt(2.0)/3.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs33, +1.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs38, -2.0/sqrt(3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs88, +1.0/3.0));


    // DSs = (2./3.)*DnS00 - 2*sqrt(2/3)*DnS03 + (2*sqrt(2)/3)*DnS08 + DnS33 - (2/sqrt(3))*DnS38 + (1/3)*DnS88 ,  where DnS is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSs = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    gsl_complex DnSs00 = DnSs.getValue(0, 0);
    gsl_complex DnSs03 = DnSs.getValue(0, 1);
    gsl_complex DnSs08 = DnSs.getValue(0, 2);
    gsl_complex DnSs33 = DnSs.getValue(1, 1);
    gsl_complex DnSs38 = DnSs.getValue(1, 2);
    gsl_complex DnSs88 = DnSs.getValue(2, 2);

    gsl_complex DSs;
    DSs = gsl_complex_mul_real(DnSs00, +2.0/3.0);
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs03, -2.0*sqrt(2.0/3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs08, +2.0*sqrt(2.0)/3.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs33, +1.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs38, -2.0/sqrt(3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs88, +1.0/3.0));


    // DPt = (2./3.)*DnP00 - 2*sqrt(2/3)*DnP03 + (2*sqrt(2)/3)*DnP08 + DnP33 - (2/sqrt(3))*DnP38 + (1/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    gsl_complex DnPt33 = DnPt.getValue(1, 1);
    gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt03, -2.0*sqrt(2.0/3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, +2.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt33, +1.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt38, -2.0/sqrt(3.0)));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +1.0/3.0));

    
    // DSt = (2./3.)*DnS00 - 2*sqrt(2/3)*DnS03 + (2*sqrt(2)/3)*DnS08 + DnS33 - (2/sqrt(3))*DnS38 + (1/3)*DnS88 ,  where DnP is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    gsl_complex DnSt33 = DnSt.getValue(1, 1);
    gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt03, -2.0*sqrt(2.0/3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, +2.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt33, +1.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt38, -2.0/sqrt(3.0)));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, +1.0/3.0));


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//DDBar->SSBar differential cross section
double dsigmadtDDBarSSBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassD;
    double m2 = effMassD;
    double m3 = effMassS;
    double m4 = effMassS;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = (2./3.)*DnP00 - sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 + (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPs = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision); 
    gsl_complex DnPs00 = DnPs.getValue(0, 0);
    gsl_complex DnPs03 = DnPs.getValue(0, 1);
    gsl_complex DnPs08 = DnPs.getValue(0, 2);
    //gsl_complex DnPs33 = DnPs.getValue(1, 1);
    gsl_complex DnPs38 = DnPs.getValue(1, 2);
    gsl_complex DnPs88 = DnPs.getValue(2, 2);

    gsl_complex DPs;
    DPs = gsl_complex_mul_real(DnPs00, +2.0/3.0);
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs03, -sqrt(2.0/3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs08, -sqrt(2.0)/3.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs38, +2.0/sqrt(3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs88, -2.0/3.0));


    // DSs = (2./3.)*DnS00 - sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 + (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSs = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    gsl_complex DnSs00 = DnSs.getValue(0, 0);
    gsl_complex DnSs03 = DnSs.getValue(0, 1);
    gsl_complex DnSs08 = DnSs.getValue(0, 2);
    //gsl_complex DnSs33 = DnSs.getValue(1, 1);
    gsl_complex DnSs38 = DnSs.getValue(1, 2);
    gsl_complex DnSs88 = DnSs.getValue(2, 2);

    gsl_complex DSs;
    DSs = gsl_complex_mul_real(DnSs00, +2.0/3.0);
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs03, -sqrt(2.0/3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs08, -sqrt(2.0)/3.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs38, +2.0/sqrt(3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs88, -2.0/3.0));


    // DPt = 2*neutralKaonPropagator
    gsl_complex DPt;
    DPt = neutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DPt = gsl_complex_mul_real(DPt, 2);


    // DSt = 2*neutralSigmaKaonPlusPropagator
    gsl_complex DSt;
    DSt = neutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DSt = gsl_complex_mul_real(DSt, 2);


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//SSBar->UUBar differential cross section
double dsigmadtSSBarUUBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassS;
    double m2 = effMassS;
    double m3 = effMassU;
    double m4 = effMassU;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = (2./3.)*DnP00 + sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 - (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPs = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision); 
    gsl_complex DnPs00 = DnPs.getValue(0, 0);
    gsl_complex DnPs03 = DnPs.getValue(0, 1);
    gsl_complex DnPs08 = DnPs.getValue(0, 2);
    //gsl_complex DnPs33 = DnPs.getValue(1, 1);
    gsl_complex DnPs38 = DnPs.getValue(1, 2);
    gsl_complex DnPs88 = DnPs.getValue(2, 2);

    gsl_complex DPs;
    DPs = gsl_complex_mul_real(DnPs00, +2.0/3.0);
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs03, +sqrt(2.0/3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs08, -sqrt(2.0)/3.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs38, -2.0/sqrt(3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs88, -2.0/3.0));
    

    // DSs = (2./3.)*DnS00 + sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 - (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSs = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    gsl_complex DnSs00 = DnSs.getValue(0, 0);
    gsl_complex DnSs03 = DnSs.getValue(0, 1);
    gsl_complex DnSs08 = DnSs.getValue(0, 2);
    //gsl_complex DnSs33 = DnSs.getValue(1, 1);
    gsl_complex DnSs38 = DnSs.getValue(1, 2);
    gsl_complex DnSs88 = DnSs.getValue(2, 2);

    gsl_complex DSs;
    DSs = gsl_complex_mul_real(DnSs00, +2.0/3.0);
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs03, +sqrt(2.0/3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs08, -sqrt(2.0)/3.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs38, -2.0/sqrt(3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs88, -2.0/3.0));


    // DPt = 2*kaonMinusPropagator
    gsl_complex DPt;
    DPt = kaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DPt = gsl_complex_mul_real(DPt, 2);


    // DSt = 2*sigmaKaonMinusPropagator
    gsl_complex DSt;
    DSt = sigmaKaonMinusPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DSt = gsl_complex_mul_real(DSt, 2);


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//SSBar->DDBar differential cross section
double dsigmadtSSBarDDBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassS;
    double m2 = effMassS;
    double m3 = effMassD;
    double m4 = effMassD;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = (2./3.)*DnP00 - sqrt(2/3)*DnP03 - (sqrt(2)/3)*DnP08 + (2/sqrt(3))*DnP38 - (2/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPs = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision); 
    gsl_complex DnPs00 = DnPs.getValue(0, 0);
    gsl_complex DnPs03 = DnPs.getValue(0, 1);
    gsl_complex DnPs08 = DnPs.getValue(0, 2);
    //gsl_complex DnPs33 = DnPs.getValue(1, 1);
    gsl_complex DnPs38 = DnPs.getValue(1, 2);
    gsl_complex DnPs88 = DnPs.getValue(2, 2);

    gsl_complex DPs;
    DPs = gsl_complex_mul_real(DnPs00, +2.0/3.0);
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs03, -sqrt(2.0/3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs08, -sqrt(2.0)/3.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs38, +2.0/sqrt(3.0)));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs88, -2.0/3.0));


    // DSs = (2./3.)*DnS00 - sqrt(2/3)*DnS03 - (sqrt(2)/3)*DnS08 + (2/sqrt(3))*DnS38 - (2/3)*DnS88 ,  where DnS is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnSs = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    gsl_complex DnSs00 = DnSs.getValue(0, 0);
    gsl_complex DnSs03 = DnSs.getValue(0, 1);
    gsl_complex DnSs08 = DnSs.getValue(0, 2);
    //gsl_complex DnSs33 = DnSs.getValue(1, 1);
    gsl_complex DnSs38 = DnSs.getValue(1, 2);
    gsl_complex DnSs88 = DnSs.getValue(2, 2);

    gsl_complex DSs;
    DSs = gsl_complex_mul_real(DnSs00, +2.0/3.0);
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs03, -sqrt(2.0/3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs08, -sqrt(2.0)/3.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs38, +2.0/sqrt(3.0)));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs88, -2.0/3.0));


    // DPt = 2*antiNeutralKaonPropagator
    gsl_complex DPt;
    DPt = antiNeutralKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DPt = gsl_complex_mul_real(DPt, 2);


    // DSt = 2*antiNeutralSigmaKaonPropagator
    gsl_complex DSt;
    DSt = antiNeutralSigmaKaonPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    DSt = gsl_complex_mul_real(DSt, 2);


    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
//SSBar->SSBar differential cross section
double dsigmadtSSBarSSBar(SU3NJL3DCutoffParameters parametersNJL, double T, 
                          double effChemPotU, double effChemPotD, double effChemPotS, 
                          double effMassU, double effMassD, double effMassS, 
                          double s, double t, double integralPrecision)
{   
    //set incoming and outgoing masses
    double m1 = effMassS;
    double m2 = effMassS;
    double m3 = effMassS;
    double m4 = effMassS;

    //calculate propagators that take part in this process
    double sK0 = sChannelZeroMomentum(s);
    double sKVec = sChannelThreeMomentum();
    double tK0 = tChannelZeroMomentum(s, m1, m2, m3, m4);
    double tKVec = tChannelThreeMomentum(s, t, m1, m2, m3, m4);

    // DPs = (2./3.)*DnP00 - (4*sqrt(2)/3)*DnP08 + (4/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPs = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision); 
    gsl_complex DnPs00 = DnPs.getValue(0, 0);
    //gsl_complex DnPs03 = DnPs.getValue(0, 1);
    gsl_complex DnPs08 = DnPs.getValue(0, 2);
    //gsl_complex DnPs33 = DnPs.getValue(1, 1);
    //gsl_complex DnPs38 = DnPs.getValue(1, 2);
    gsl_complex DnPs88 = DnPs.getValue(2, 2);

    gsl_complex DPs;
    DPs = gsl_complex_mul_real(DnPs00, +2.0/3.0);
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs08, -4.0*sqrt(2.0)/3.0));
    DPs = gsl_complex_add(DPs, gsl_complex_mul_real(DnPs88, +4.0/3.0));


    // DSs = (2./3.)*DnS00 - (4*sqrt(2)/3)*DnS08 + (4/3)*DnS88 ,  where DnS is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSs = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, sK0, sKVec, integralPrecision);
    gsl_complex DnSs00 = DnSs.getValue(0, 0);
    //gsl_complex DnSs03 = DnSs.getValue(0, 1);
    gsl_complex DnSs08 = DnSs.getValue(0, 2);
    //gsl_complex DnSs33 = DnSs.getValue(1, 1);
    //gsl_complex DnSs38 = DnSs.getValue(1, 2);
    gsl_complex DnSs88 = DnSs.getValue(2, 2);

    gsl_complex DSs;
    DSs = gsl_complex_mul_real(DnSs00, +2.0/3.0);
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs08, -4.0*sqrt(2.0)/3.0));
    DSs = gsl_complex_add(DSs, gsl_complex_mul_real(DnSs88, +4.0/3.0));


    // DPt = (2./3.)*DnP00 - (4*sqrt(2)/3)*DnP08 + (4/3)*DnP88 ,  where DnP is the neutral pseudoscalar propagator matrix
    ComplexSquareMatrixGSL DnPt = 
    neutral038PseudoscalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision); 
    gsl_complex DnPt00 = DnPt.getValue(0, 0);
    //gsl_complex DnPt03 = DnPt.getValue(0, 1);
    gsl_complex DnPt08 = DnPt.getValue(0, 2);
    //gsl_complex DnPt33 = DnPt.getValue(1, 1);
    //gsl_complex DnPt38 = DnPt.getValue(1, 2);
    gsl_complex DnPt88 = DnPt.getValue(2, 2);

    gsl_complex DPt;
    DPt = gsl_complex_mul_real(DnPt00, +2.0/3.0);
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt08, -4.0*sqrt(2.0)/3.0));
    DPt = gsl_complex_add(DPt, gsl_complex_mul_real(DnPt88, +4.0/3.0));

    
    // DSt = (2./3.)*DnS00 - (4*sqrt(2)/3)*DnS08 + (4/3)*DnS88 ,  where DnP is the neutral scalar propagator matrix
    ComplexSquareMatrixGSL DnSt = 
    neutral038ScalarsPropagator(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, tK0, tKVec, integralPrecision);
    gsl_complex DnSt00 = DnSt.getValue(0, 0);
    //gsl_complex DnSt03 = DnSt.getValue(0, 1);
    gsl_complex DnSt08 = DnSt.getValue(0, 2);
    //gsl_complex DnSt33 = DnSt.getValue(1, 1);
    //gsl_complex DnSt38 = DnSt.getValue(1, 2);
    gsl_complex DnSt88 = DnSt.getValue(2, 2);

    gsl_complex DSt;
    DSt = gsl_complex_mul_real(DnSt00, +2.0/3.0);
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt08, -4.0*sqrt(2.0)/3.0));
    DSt = gsl_complex_add(DSt, gsl_complex_mul_real(DnSt88, +4.0/3.0));

    //calculate the differential cross section using the above propagators
    double Nc = parametersNJL.getNumberOfColours();
    double dsigmadt = dsigmadtNJLQuarkAntiquarkScattering(Nc, s, t, m1, m2, m3, m4, DPs, DSs, DPt, DSt);

    return dsigmadt;
}


////////////////////////////////////////////////////////////////////////////////////////
double differentialCrossSectionProcess12To34(SU3NJL3DCutoffParameters parametersNJL, double T, 
                                             double effChemPotU, double effChemPotD, double effChemPotS, 
                                             double effMassU, double effMassD, double effMassS, 
                                             double s, double t, double integralPrecision, scatteringProcess process)
{
    double dsigmadtProcess = 0.0;

    if ( process==UDUD )
    {
        dsigmadtProcess = 
        dsigmadtUDUD(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);
    }
    else if( process==DUDU )
    {
        dsigmadtProcess = 
        dsigmadtDUDU(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);
    }
    else if( process==USUS )
    {
        dsigmadtProcess = 
        dsigmadtUSUS(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);     
    }
    else if( process==SUSU )
    {
        dsigmadtProcess = 
        dsigmadtSUSU(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);         
    }
    else if( process==DSDS )
    {
        dsigmadtProcess = 
        dsigmadtDSDS(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);         
    }
    else if( process==SDSD )
    {
        dsigmadtProcess = 
        dsigmadtSDSD(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);           
    }
    else if( process==UUUU )
    {
        dsigmadtProcess = 
        dsigmadtUUUU(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==DDDD )
    {
        dsigmadtProcess = 
        dsigmadtDDDD(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);       
    }
    else if( process==SSSS )
    {
        dsigmadtProcess = 
        dsigmadtSSSS(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);   
    }
    else if( process==UDBarUDBar )
    {
        dsigmadtProcess = 
        dsigmadtUDBarUDBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);     
    }
    else if( process==USBarUSBar )
    {
        dsigmadtProcess = 
        dsigmadtUSBarUSBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);           
    }
    else if( process==DSBarDSBar )
    {
        dsigmadtProcess = 
        dsigmadtDSBarDSBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);           
    }
    else if( process==DUBarDUBar )
    {
        dsigmadtProcess = 
        dsigmadtDUBarDUBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);              
    }
    else if( process==SUBarSUBar )
    {
        dsigmadtProcess = 
        dsigmadtSUBarSUBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);              
    }
    else if( process==SDBarSDBar )
    {
        dsigmadtProcess = 
        dsigmadtSDBarSDBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);      
    }
    else if( process==UUBarUUBar )
    {
        dsigmadtProcess = 
        dsigmadtUUBarUUBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==UUBarDDBar )
    {
        dsigmadtProcess = 
        dsigmadtUUBarDDBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==UUBarSSBar )
    {
        dsigmadtProcess = 
        dsigmadtUUBarSSBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==DDBarUUBar )
    {
        dsigmadtProcess = 
        dsigmadtDDBarUUBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);         
    }
    else if( process==DDBarDDBar )
    {
        dsigmadtProcess = 
        dsigmadtDDBarDDBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);         
    }
    else if( process==DDBarSSBar )
    {
        dsigmadtProcess = 
        dsigmadtDDBarSSBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);    
    }
    else if( process==SSBarUUBar )
    {
        dsigmadtProcess = 
        dsigmadtSSBarUUBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==SSBarDDBar )
    {
        dsigmadtProcess = 
        dsigmadtSSBarDDBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==SSBarSSBar )
    {
        dsigmadtProcess = 
        dsigmadtSSBarSSBar(parametersNJL, T, effChemPotU, effChemPotD, effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==UBarUBarUBarUBar )
    {
        dsigmadtProcess = 
        dsigmadtUUUU(parametersNJL, T, -effChemPotU, -effChemPotD, -effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==UBarDBarUBarDBar )
    {
        dsigmadtProcess = 
        dsigmadtUDUD(parametersNJL, T, -effChemPotU, -effChemPotD, -effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==UBarSBarUBarSBar )
    {
        dsigmadtProcess = 
        dsigmadtUSUS(parametersNJL, T, -effChemPotU, -effChemPotD, -effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else if( process==SBarSBarSBarSBar )
    {
        dsigmadtProcess = 
        dsigmadtSSSS(parametersNJL, T, -effChemPotU, -effChemPotD, -effChemPotS, effMassU, effMassD, effMassS, s, t, integralPrecision);            
    }
    else
    {
        cout << "Differential cross section for the provided process is not defined in differentialCrossSectionProcess12To34! Aborting!\n";
        abort();
    }

    return dsigmadtProcess;
}


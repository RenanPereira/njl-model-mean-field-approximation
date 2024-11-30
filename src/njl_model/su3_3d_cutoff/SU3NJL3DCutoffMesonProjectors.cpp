#include <cmath>
#include <iostream>
#include "OneFermionLineIntegral.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffMesonProjectors.h"

using namespace std;



////////////////////////////////////////////////////////////////////////////////////////
//SU3 NJL model Pseudoscalar Meson Projectors


double SU3NJL3DCutoffMesonProjector::pseudoscalar00()
{
	double P00 = ( 1.0/2.0 )*( + gS 
							   - ( 2.0/3.0 )*kappa*( sigmaU + sigmaD + sigmaS )
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + ( 2.0/3.0 )*g2*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) ) );
	return P00;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar08()
{
	double P08 = ( 1.0/2.0 )*( + ( kappa/(3.0*sqrt(2)) )*( sigmaU + sigmaD - 2*sigmaS )
							   + ( g2*sqrt(2.0)/3.0 )*( pow(sigmaU,2) + pow(sigmaD,2) - 2*pow(sigmaS,2) ) );

	return P08;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar03()
{
	double P03 = ( 1.0/2.0 )*( + ( kappa/( sqrt(6.0) ) )*( sigmaU - sigmaD )
							   + ( g2*sqrt(2.0/3.0) )*( pow(sigmaU,2) - pow(sigmaD,2)  ) );

	return P03;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar33()
{
	double P33 = ( 1.0/2.0 )*( + gS
							   + kappa*sigmaS
			     			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
			     			   + g2*( pow(sigmaU,2) + pow(sigmaD,2) ) );

	return P33;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar38()
{
	double P38 = ( 1.0/2.0 )*( - ( kappa/( sqrt(3.0) ) )*( sigmaU - sigmaD )
							   + ( g2/sqrt(3.0) )*( pow(sigmaU,2) - pow(sigmaD,2)  ) );

	return P38;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar88()
{
	double P88 = ( 1.0/2.0 )*( + gS
							   + ( 1.0/3.0 )*kappa*( 2.0*sigmaU + 2.0*sigmaD - sigmaS )
			     	           + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
			     			   + ( 1.0/3.0 )*g2*( pow(sigmaU,2) + pow(sigmaD,2) + 4.0*pow(sigmaS,2) ) );

	return P88;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar11()
{
	double P11 = ( 1.0/2.0 )*( + gS
						       + kappa*sigmaS
			     			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
			     			   + 2.0*g2*( pow(sigmaU,2) + pow(sigmaD,2) - sigmaU*sigmaD ) );

	return P11;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar22()
{
	double P22 = ( 1.0/2.0 )*( + gS
						       + kappa*sigmaS
			     			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
			     			   + 2.0*g2*( pow(sigmaU,2) + pow(sigmaD,2) - sigmaU*sigmaD ) );

	return P22;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar44()
{
	double P44 = ( 1.0/2.0 )*( + gS
							   + kappa*sigmaD
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaU,2) + pow(sigmaS,2) - sigmaU*sigmaS ) );

	return P44;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar55()
{
	double P55 = ( 1.0/2.0 )*( + gS
							   + kappa*sigmaD
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaU,2) + pow(sigmaS,2) - sigmaU*sigmaS ) );

	return P55;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar66()
{
	double P66 = ( 1.0/2.0 )*( + gS
							   + kappa*sigmaU
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaD,2) + pow(sigmaS,2) - sigmaD*sigmaS ) );

	return P66;
}


double SU3NJL3DCutoffMesonProjector::pseudoscalar77()
{
	double P77 = ( 1.0/2.0 )*( + gS
							   + kappa*sigmaU
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaD,2) + pow(sigmaS,2) - sigmaD*sigmaS ) );

	return P77;
}


////////////////////////////////////////////////////////////////////////////////////////
//SU3 NJL model Scalar Meson Projectors


double SU3NJL3DCutoffMesonProjector::scalar00()
{
	double S00 = ( 1.0/2.0 )*( + gS
							   + ( 2.0/3.0 )*kappa*( sigmaU + sigmaD + sigmaS )
				 			   + ( 2.0/3.0 )*g1*( 5.0*pow(sigmaU,2) + 5.0*pow(sigmaD,2) + 5.0*pow(sigmaS,2) + 4.0*sigmaU*sigmaS + 4.0*sigmaS*sigmaD + 4.0*sigmaU*sigmaD )
				 			   + 2.0*g2*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) ) );

	return S00;
}


double SU3NJL3DCutoffMesonProjector::scalar08()
{
	double S08 = ( 1.0/2.0 )*( + sqrt(2.0)*( -1.0/6.0 )*kappa*( sigmaU + sigmaD - 2*sigmaS )
			   				   + sqrt(2.0)*( 2.0/3.0 )*g1*( pow(sigmaU,2) + pow(sigmaD,2) - 2.0*pow(sigmaS,2) + 2.0*sigmaU*sigmaD - sigmaU*sigmaS - sigmaD*sigmaS )
			   				   + sqrt(2.0)*g2*( pow(sigmaU,2) + pow(sigmaD,2) - 2.0*pow(sigmaS,2) ) );

	return S08;
}


double SU3NJL3DCutoffMesonProjector::scalar03()
{
	double S03 = ( 1.0/2.0 )*( - ( kappa/sqrt(6.0) )*( sigmaU - sigmaD )
				 			   + 2.0*sqrt(2.0/3.0)*g1*( sigmaU - sigmaD )*( sigmaU + sigmaD + sigmaS )
				 			   + sqrt(6.0)*g2*( pow(sigmaU,2) - pow(sigmaD,2)  ) );

	return S03;
}


double SU3NJL3DCutoffMesonProjector::scalar38()
{
	double S38 = ( 1.0/2.0 )*( + ( kappa/sqrt(3.0) )*( sigmaU - sigmaD )
				 			   + ( 2.0/sqrt(3.0) )*g1*( sigmaU - sigmaD )*( sigmaU + sigmaD - 2.0*sigmaS )
				 			   + sqrt(3.0)*g2*( pow(sigmaU,2) - pow(sigmaD,2)  ) );

	return S38;
}


double SU3NJL3DCutoffMesonProjector::scalar33()
{
	double S33 = ( 1.0/2.0 )*( + gS
							   - kappa*sigmaS
				 			   + 4.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + (1.0/2.0)*pow(sigmaS,2) - sigmaU*sigmaD )
				 			   + 3.0*g2*( pow(sigmaU,2) + pow(sigmaD,2) ) );

	return S33;
}


double SU3NJL3DCutoffMesonProjector::scalar88()
{
	double S88 = ( 1.0/2.0 )*( + gS
							   - ( 1.0/3.0 )*kappa*( 2.0*sigmaU + 2.0*sigmaD - sigmaS )
			     			   + ( 2.0/3.0 )*g1*( 4.0*pow(sigmaU,2) + 4.0*pow(sigmaD,2) + 7.0*pow(sigmaS,2) + 2.0*sigmaU*sigmaD - 4.0*sigmaU*sigmaS - 4.0*sigmaD*sigmaS )
				 			   + g2*( pow(sigmaU,2) + pow(sigmaD,2) + 4*pow(sigmaS,2) ) );

	return S88;
}


double SU3NJL3DCutoffMesonProjector::scalar11()
{
	double S11 = ( 1.0/2.0 )*( + gS
							   - kappa*sigmaS
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaU,2) + pow(sigmaD,2) + sigmaU*sigmaD ) );

	return S11;
}


double SU3NJL3DCutoffMesonProjector::scalar22()
{
	double S22 = ( 1.0/2.0 )*( + gS
							   - kappa*sigmaS
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaU,2) + pow(sigmaD,2) + sigmaU*sigmaD ) );

	return S22;
}


double SU3NJL3DCutoffMesonProjector::scalar44()
{
	double S44 = ( 1.0/2.0 )*( + gS
							   - kappa*sigmaD
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaU,2) + pow(sigmaS,2) + sigmaU*sigmaS ) );

	return S44;
}


double SU3NJL3DCutoffMesonProjector::scalar55()
{
	double S55 = ( 1.0/2.0 )*( + gS
							   - kappa*sigmaD
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaU,2) + pow(sigmaS,2) + sigmaU*sigmaS ) );

	return S55;
}


double SU3NJL3DCutoffMesonProjector::scalar66()
{
	double S66 = ( 1.0/2.0 )*( + gS
							   - kappa*sigmaU
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaD,2) + pow(sigmaS,2) + sigmaD*sigmaS ) );

	return S66;
}


double SU3NJL3DCutoffMesonProjector::scalar77()
{
	double S77 = ( 1.0/2.0 )*( + gS
							   - kappa*sigmaU
				 			   + 2.0*g1*( pow(sigmaU,2) + pow(sigmaD,2) + pow(sigmaS,2) )
				 			   + 2.0*g2*( pow(sigmaD,2) + pow(sigmaS,2) + sigmaD*sigmaS ) );

	return S77;
}




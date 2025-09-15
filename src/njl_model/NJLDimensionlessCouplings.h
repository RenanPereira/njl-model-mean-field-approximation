#ifndef NJLDIMENSIONLESSCOUPLINGS_H
#define NJLDIMENSIONLESSCOUPLINGS_H

#include <string>


namespace NJLDimensionlessCouplings 
{
    inline const std::string FOUR_QUARK_SP_COUPLING = "fourQuarkSPCouplingCutoff2";//g
	inline const std::string FOUR_QUARK_VP_COUPLING = "fourQuarkVPCouplingCutoff";//gOmega
	inline const std::string FOUR_QUARK_VIPI_COUPLING = "fourQuarkVIPICouplingCutoff";//gRho

    inline const std::string DETERMINANT_COUPLING = "determinantCouplingCutoff5";//kappa
	
    inline const std::string EIGHT_QUARK_SP_OZI_COUPLING = "eightQuarkSPOziViolatingCouplingCutoff8";//g1
    inline const std::string EIGHT_QUARK_SP_NON_OZI_COUPLING = "eightQuarkSPNonOziViolatingCouplingCutoff8";//g2
	inline const std::string EIGHT_QUARK_VP_COUPLING = "eightQuarkVPCoupling";//gOmegaOmega
	inline const std::string EIGHT_QUARK_VIPI_COUPLING = "eightQuarkVIPICoupling";//gRhoRho
	inline const std::string EIGHT_QUARK_VPVIPI_COUPLING = "eightQuarkVPVIPICoupling";//gOmegaRho
	inline const std::string EIGHT_QUARK_SPVP_COUPLING = "eightQuarkSPVPCoupling";//gSigmaOmega
	inline const std::string EIGHT_QUARK_SPVIPI_COUPLING = "eightQuarkSPVIPICoupling";//gSigmaRho

	inline const std::string TWELVE_QUARK_VP_COUPLING = "twelveQuarkVPCoupling";//gOmega3
	
	inline const std::string SIXTEEN_QUARK_VP_COUPLING = "sixteenQuarkVPCoupling";//gOmega4
}


#endif
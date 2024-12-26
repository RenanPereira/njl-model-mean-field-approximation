#ifndef NJLDIMENSIONLESSCOUPLINGS_H
#define NJLDIMENSIONLESSCOUPLINGS_H

#include <string>


namespace NJLDimensionlessCouplings {
    const std::string FOUR_QUARK_SP_COUPLING = "fourQuarkSPCouplingCutoff2";//g
	const std::string FOUR_QUARK_VP_COUPLING = "fourQuarkVPCouplingCutoff";//gOmega
	const std::string FOUR_QUARK_VIPI_COUPLING = "fourQuarkVIPICouplingCutoff";//gRho

    const std::string DETERMINANT_COUPLING = "determinantCouplingCutoff5";//kappa
	
    const std::string EIGHT_QUARK_SP_OZI_COUPLING = "eightQuarkSPOziViolatingCouplingCutoff8";//g1
    const std::string EIGHT_QUARK_SP_NON_OZI_COUPLING = "eightQuarkSPNonOziViolatingCouplingCutoff8";//g2
	const std::string EIGHT_QUARK_VP_COUPLING = "eightQuarkVPCoupling";//gOmegaOmega
	const std::string EIGHT_QUARK_VIPI_COUPLING = "eightQuarkVIPICoupling";//gRhoRho
	const std::string EIGHT_QUARK_VPVIPI_COUPLING = "eightQuarkVPVIPICoupling";//gOmegaRho
	const std::string EIGHT_QUARK_SPVP_COUPLING = "eightQuarkSPVPCoupling";//gSigmaOmega
	const std::string EIGHT_QUARK_SPVIPI_COUPLING = "eightQuarkSPVIPICoupling";//gSigmaRho

	const std::string TWELVE_QUARK_VP_COUPLING = "twelveQuarkVPCoupling";//gOmega3
	
	const std::string SIXTEEN_QUARK_VP_COUPLING = "sixteenQuarkVPCoupling";//gOmega4
}


#endif
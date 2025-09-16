#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H


namespace PhysicalConstants
{
    inline constexpr double hc_MeVfm = 197.3269804;
    inline constexpr double hc_GeVfm = hc_MeVfm/1000.0;

    inline constexpr double electronMass_MeV = 0.511;
    inline constexpr double electronMass_GeV = electronMass_MeV/1000.0;

    // 1 bn = 2568 GeV^{−2} 
    // 1 GeV^{−2} = 3.894x10^{−4} bn
    inline constexpr double barnToInverseGeVSquared = 2568.0;
    inline constexpr double inverseGeVSquaredToMiliBarn = 0.3894;
    inline constexpr double inverseGeVSquaredToBarn = inverseGeVSquaredToMiliBarn*1E-3;
}


#endif

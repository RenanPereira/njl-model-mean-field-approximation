#ifndef NJL_REGULARIZATION_SCHEMES_H
#define NJL_REGULARIZATION_SCHEMES_H

#include <map>
#include <string>


// Enum representing various 3D cutoff regularization schemes for the NJL model
enum class NJL3DCutoffRegularizationScheme 
{
    CUTOFF_EVERYWHERE,                    // Apply cutoff everywhere
    CUTOFF_EVERYWHERE_WITH_CTMU,          // Apply cutoff everywhere with CTmu
    CUTOFF_ON_DIVERGENT_INTEGRALS_ONLY    // Apply cutoff only on divergent integrals
};

inline const std::map<NJL3DCutoffRegularizationScheme, std::string> NJL3DCutoffRegularizationSchemeMap = 
{
    {NJL3DCutoffRegularizationScheme::CUTOFF_EVERYWHERE, "CUTOFF_EVERYWHERE"},
    {NJL3DCutoffRegularizationScheme::CUTOFF_EVERYWHERE_WITH_CTMU, "CUTOFF_EVERYWHERE_WITH_CTMU"},
    {NJL3DCutoffRegularizationScheme::CUTOFF_ON_DIVERGENT_INTEGRALS_ONLY, "CUTOFF_ON_DIVERGENT_INTEGRALS_ONLY"}
};

std::string toString(NJL3DCutoffRegularizationScheme );

NJL3DCutoffRegularizationScheme stringToNJL3DCutoffRegularizationScheme(const std::string& );

bool isValidNJL3DCutoffRegularizationScheme(const std::string& );

#endif

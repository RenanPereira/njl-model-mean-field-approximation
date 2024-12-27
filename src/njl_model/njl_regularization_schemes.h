#ifndef NJL_REGULARIZATION_SCHEMES_H
#define NJL_REGULARIZATION_SCHEMES_H

#include <map>


// Enum representing various 3D cutoff regularization schemes for the NJL model
enum NJL3DCutoffRegularizationScheme 
{
    cutoffEverywhere,                         // Apply cutoff everywhere
    cutoffEverywhereWithCTmu,                 // Apply cutoff everywhere with CTmu
    cutoffOnDivergentIntegralsOnly            // Apply cutoff only on divergent integrals
};

static const std::map<NJL3DCutoffRegularizationScheme, std::string> NJL3DCutoffRegularizationSchemeMap = 
{
    {cutoffEverywhere, "cutoffEverywhere"},
    {cutoffEverywhereWithCTmu, "cutoffEverywhereWithCTmu"},
    {cutoffOnDivergentIntegralsOnly, "cutoffOnDivergentIntegralsOnly"}
};

std::string toStringNJL3DCutoffRegularizationScheme(NJL3DCutoffRegularizationScheme );

NJL3DCutoffRegularizationScheme stringToNJL3DCutoffRegularizationScheme(const std::string& );

bool isValidNJL3DCutoffRegularizationScheme(const std::string& );

#endif

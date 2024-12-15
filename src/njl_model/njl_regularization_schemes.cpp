#include <iostream>
#include <string>
#include "njl_model/njl_regularization_schemes.h"

using namespace std;

string toStringNJL3DCutoffRegularizationScheme(NJL3DCutoffRegularizationScheme scheme) 
{
    switch (scheme) 
    {
        case cutoffEverywhere:
            return "cutoffEverywhere";
        case cutoffEverywhereWithCTmu:
            return "cutoffEverywhereWithCTmu";
        case cutoffOnDivergentIntegralsOnly:
            return "cutoffOnDivergentIntegralsOnly";
        default:
            cout << "Invalid NJL3DCutoffRegularizationScheme value, returning Unknown" << endl;
            return "Unknown";
    }
}


NJL3DCutoffRegularizationScheme stringToNJL3DCutoffRegularizationScheme(const string& schemeStr) 
{
    if (schemeStr == "cutoffEverywhere") 
    {
        return cutoffEverywhere;
    } 
    else if (schemeStr == "cutoffEverywhereWithCTmu") 
    {
        return cutoffEverywhereWithCTmu;
    } 
    else if (schemeStr == "cutoffOnDivergentIntegralsOnly") 
    {
        return cutoffOnDivergentIntegralsOnly;
    } 
    else 
    {
        cout << "Invalid NJL3DCutoffRegularizationScheme string: " + schemeStr + ". Aborting!\n";
        abort();
    }
}
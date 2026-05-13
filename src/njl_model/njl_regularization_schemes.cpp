#include <iostream>
#include <string>
#include "njl_model/njl_regularization_schemes.h"


std::string toString(NJL3DCutoffRegularizationScheme scheme) 
{
    // Check if the method exists in the map using count
    if (NJL3DCutoffRegularizationSchemeMap.count(scheme))
    {
        return NJL3DCutoffRegularizationSchemeMap.at(scheme);
    } 
    else 
    {
        std::cout << "Error: NJL3DCutoffRegularizationScheme not found in map! Returning UNKNOWN." << std::endl;
        return "UNKNOWN";
    }
}


NJL3DCutoffRegularizationScheme stringToNJL3DCutoffRegularizationScheme(const std::string& schemeString) 
{
    // Iterate over the map with explicit type
    for (std::map<NJL3DCutoffRegularizationScheme, std::string>::const_iterator it = NJL3DCutoffRegularizationSchemeMap.begin(); it != NJL3DCutoffRegularizationSchemeMap.end(); ++it) 
    {
        if (it->second == schemeString) 
        {
            return it->first;
        }
    }

    std::cout << "Invalid NJL3DCutoffRegularizationScheme string: " + schemeString + ". Aborting!\n";
    abort();
}


bool isValidNJL3DCutoffRegularizationScheme(const std::string& schemeString)
{
    bool isRegularizationSchemeValid = false;
    // Iterate over the map with explicit type
    for (std::map<NJL3DCutoffRegularizationScheme, std::string>::const_iterator it = NJL3DCutoffRegularizationSchemeMap.begin(); it != NJL3DCutoffRegularizationSchemeMap.end(); ++it) 
    {
        if (it->second == schemeString) 
        {
            isRegularizationSchemeValid = true;
            break;
        }
    }

    if( isRegularizationSchemeValid==false )
    {
        std::cout << "The value " + schemeString + " is not a NJL3DCutoffRegularizationScheme!\n";
    }
    
    return isRegularizationSchemeValid;
}

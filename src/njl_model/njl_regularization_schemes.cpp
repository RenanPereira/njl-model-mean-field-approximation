#include <iostream>
#include <string>
#include "njl_model/njl_regularization_schemes.h"

using namespace std;


string toStringNJL3DCutoffRegularizationScheme(NJL3DCutoffRegularizationScheme scheme) 
{
    // Check if the method exists in the map using count
    if (NJL3DCutoffRegularizationSchemeMap.count(scheme))
    {
        return NJL3DCutoffRegularizationSchemeMap.at(scheme);
    } 
    else 
    {
        cout << "Error: NJL3DCutoffRegularizationScheme not found in map! Returning UNKNOWN." << endl;
        return "UNKNOWN";
    }
}


NJL3DCutoffRegularizationScheme stringToNJL3DCutoffRegularizationScheme(const string& schemeString) 
{
    // Iterate over the map with explicit type
    for (map<NJL3DCutoffRegularizationScheme, string>::const_iterator it = NJL3DCutoffRegularizationSchemeMap.begin(); it != NJL3DCutoffRegularizationSchemeMap.end(); ++it) 
    {
        if (it->second == schemeString) 
        {
            return it->first;
        }
    }

    cout << "Invalid NJL3DCutoffRegularizationScheme string: " + schemeString + ". Aborting!\n";
    abort();
}


bool isValidNJL3DCutoffRegularizationScheme(const string& schemeString)
{
    bool isRegularizationSchemeValid = false;
    // Iterate over the map with explicit type
    for (map<NJL3DCutoffRegularizationScheme, string>::const_iterator it = NJL3DCutoffRegularizationSchemeMap.begin(); it != NJL3DCutoffRegularizationSchemeMap.end(); ++it) 
    {
        if (it->second == schemeString) 
        {
            isRegularizationSchemeValid = true;
            break;
        }
    }

    if( isRegularizationSchemeValid==false )
    {
        cout << "The value " + schemeString + " is not a NJL3DCutoffRegularizationScheme!\n";
    }
    
    return isRegularizationSchemeValid;
}

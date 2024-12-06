#ifndef NJL_REGULARIZATION_SCHEMES_H
#define NJL_REGULARIZATION_SCHEMES_H


// Enum representing various 3D cutoff regularization schemes for the NJL model
enum NJL3DCutoffRegularizationScheme {
    cutoffEverywhere,                // Apply cutoff everywhere
    cutoffEverywhereWithCTmu,        // Apply cutoff everywhere with CTmu
    cutoffOnDivergentIntegralsOnly   // Apply cutoff only on divergent integrals
};


#endif

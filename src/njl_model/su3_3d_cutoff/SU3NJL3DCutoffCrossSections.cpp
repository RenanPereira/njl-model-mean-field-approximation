#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include "physics_utils/distribution_functions.h"
#include "physics_utils/physical_constants.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffCrossSections.h"
#include "njl_model/su3_3d_cutoff/SU3NJL3DCutoffFixedChemPotTemp.h"


//In this file we have functions to calculate quark-quark and quark-antiquark cross sections
//quark-quark: 
//UD->UD, DU->DU, US->US, SU->SU, DS->DS, SD->SD, UU->UU, DD->DD, SS->SS
//quark-antiquark: 
//UDBar->UDBar, USBar->USBar, DSBar->DSBar,
//DUBar->DUBar, SUBar->SUBar, SDBar->SDBar,
//UUBar->UUBar, UUBar->DDBar, UUBar->SSBar,
//DDBar->UUBar, DDBar->DDBar, DDBar->SSBar,
//SSBar->UUBar, SSBar->DDBar, SSBar->SSBar
//UBarUBar->UBarUBar, UBarDBar->UBarDBar, UBarSBar->UBarSBar, SBarSBar->SBarSBar


////////////////////////////////////////////////////////////////////////////////////////
//kinematics


//momentum in the center of mass
double momentumCM(double s, double m1, double m2)
{   
    double sqrtArg = ( s - pow(m1+m2,2) )*( s - pow(m1-m2,2) );

    double pCM = sqrt( sqrtArg )/( 2.0*sqrt(s) );

    return pCM;
}


double energyOutParticle1(double s, double m1, double m2)
{
    double E1 = (1.0/2.0)*sqrt( pow(+pow(m1,2) - pow(m2,2) + s, 2)/s );

    return E1;
}


double energyOutParticle2(double s, double m1, double m2)
{
    double E2 = (1.0/2.0)*sqrt( pow(-pow(m1,2) + pow(m2,2) + s, 2)/s );

    return E2;
}


double energyOutParticle3(double s, double m3, double m4)
{
    double E3 = (1.0/2.0)*sqrt( pow(+pow(m3,2) - pow(m4,2) + s, 2)/s );

    return E3;
}


double energyOutParticle4(double s, double m3, double m4)
{
    double E4 = (1.0/2.0)*sqrt( pow(-pow(m3,2) + pow(m4,2) + s, 2)/s );

    return E4;
}


//relative velocity between two particles in the initial state (in the center of mass)
double relativeVelocityCM(double s, double m1, double m2)
{   
    double pCM = momentumCM(s, m1, m2);
    double E1 = energyOutParticle1(s, m1, m2);
    double E2 = energyOutParticle2(s, m1, m2);

    double vRel = pCM/E1 + pCM/E2;

    return vRel;
}


//relative velocity between two particles in the initial state
double relativeVelocity(double p1, double p2, double E1, double E2, double theta)
{   
    //caculate velocities from momenta and energies
    double v1 = p1/E1;
    double v2 = p2/E2;

    double vRel = sqrt( pow(v1,2) + pow(v2,2) - 2.0*v1*v2*cos(theta) - pow(v1,2)*pow(v2,2)*pow(sin(theta),2) );

    return vRel;
}


////////////////////////////////////////////////////////////////////////////////////////
//minimum and maximum values for the t channel


double tChannelMin(double s, double m1, double m2, double m3, double m4)
{   
    double sqrt12 = sqrt( pow(s + pow(m1,2) - pow(m2,2), 2)/( 4.0*s ) - pow(m1,2) );
    double sqrt34 = sqrt( pow(s + pow(m3,2) - pow(m4,2), 2)/( 4.0*s ) - pow(m3,2) );

    double tMin = pow(m1,2) + pow(m3,2) - (s + pow(m1,2) - pow(m2,2))*(s + pow(m3,2) - pow(m4,2))/( 2.0*s ) - 2.0*sqrt12*sqrt34;

    return tMin;
}


double tChannelMax(double s, double m1, double m2, double m3, double m4)
{   
    double sqrt12 = sqrt( pow(s + pow(m1,2) - pow(m2,2), 2)/( 4.0*s ) - pow(m1,2) );
    double sqrt34 = sqrt( pow(s + pow(m3,2) - pow(m4,2), 2)/( 4.0*s ) - pow(m3,2) );

    double tMax = pow(m1,2) + pow(m3,2) - (s + pow(m1,2) - pow(m2,2))*(s + pow(m3,2) - pow(m4,2))/( 2.0*s ) + 2.0*sqrt12*sqrt34;

    return tMax;
}


double centerOfMassEnergyThreshold(double m1, double m2, double m3, double m4)
{
    double m1_plus_m2_squared = pow(m1+m2, 2);
    double m3_plus_m4_squared = pow(m3+m4, 2);

    double sMin;
    if ( m1_plus_m2_squared>m3_plus_m4_squared ){ sMin = m1_plus_m2_squared; }
    else{ sMin = m3_plus_m4_squared; }

    return sMin;
}


double sMaximum(double cutoff, double m1, double m2, double m3, double m4)
{
    double sMax12 = pow( sqrt( pow(cutoff,2) + pow(m1,2) ) + sqrt( pow(cutoff,2) + pow(m2,2) ) , 2);
    double sMax34 = pow( sqrt( pow(cutoff,2) + pow(m3,2) ) + sqrt( pow(cutoff,2) + pow(m4,2) ) , 2);

    double sMax;
    if ( sMax12<sMax34 ){ sMax = sMax12; }
    else{ sMax = sMax34; }

    return sMax;
}

double sMaximumKlevansky(double cutoff, double Mu, double Md, double Ms)
{   
    double M;
    if ( Mu<=Md && Mu<=Ms ){ M = Mu; }
    else if ( Md<=Mu && Md<=Ms ){ M = Md; }
    else{ M = Ms; }

    double sMax = pow( 2.0*sqrt( pow(cutoff,2) + pow(M,2) ), 2);

    return sMax;
}


double crossSectionProcess12To34Integrand(double x, void *parameters)
{   
    CrossSectionIntegrand aux(parameters);
    SU3NJL3DCutoffParameters parametersNJL = aux.getParametersNJL();
    double T = aux.getTemperature();
    double effChemPotU = aux.getUpQuarkEffectiveChemicalPotential();
    double effChemPotD = aux.getDownQuarkEffectiveChemicalPotential();
    double effChemPotS = aux.getStrangeQuarkEffectiveChemicalPotential();
    double effMassU = aux.getUpQuarkEffectiveMass();
    double effMassD = aux.getDownQuarkEffectiveMass();
    double effMassS = aux.getStrangeQuarkEffectiveMass();
    double s = aux.getCenterOfMassEnergy();
    double integralPrecision = aux.getPropagatorIntegralPrecision();
    scatteringProcess process = aux.getProcess();
    bool largeAngleScatteringContribution = aux.getLargeAngleScatteringContribution();

    double integrand = 
    differentialCrossSectionProcess12To34(parametersNJL, T, 
                                          effChemPotU, effChemPotD, effChemPotS, 
                                          effMassU,    effMassD,    effMassS, 
                                          s, x*s, integralPrecision, process);

    if ( largeAngleScatteringContribution==true )
    {   
        double m1, m2, m3, m4;
        inOutMassesGivenScatteringProcess(effMassU, effMassD, effMassS, process, m1, m2, m3, m4);
        double sinTheta = sinScatteringAngle(s, x*s, m1, m2, m3, m4);
        integrand = integrand*pow(sinTheta, 2);
    }

    return integrand;
}


double crossSectionProcess12To34(SU3NJL3DCutoffParameters parametersNJL, double T, 
                                 double effChemPotU, double effChemPotD, double effChemPotS, 
                                 double effMassU, double effMassD, double effMassS, 
                                 double s, double propIntPrecision, scatteringProcess process, 
                                 bool largeAngleScatteringContribution, double crossSecIntPrecision)
{   
    string integralID = "crossSectionProcess" + toString(process);
    CrossSectionIntegrand aux(integralID, parametersNJL, T, 
                              effChemPotU, effChemPotD, effChemPotS, 
                              effMassU, effMassD, effMassS, 
                              s, propIntPrecision, process, largeAngleScatteringContribution);

    //set incoming and outgoing masses
    double m1, m2, m3, m4;
    inOutMassesGivenScatteringProcess(effMassU, effMassD, effMassS, process, m1, m2, m3, m4);

    //set outgoing chemical potentials
    double cP3, cP4;
    outChemicalPotentialsGivenScatteringProcess(effChemPotU, effChemPotD, effChemPotS, process, cP3, cP4);

    //calculate the energy threshold for the reaction
    double sMin = centerOfMassEnergyThreshold(m1, m2, m3, m4);
    if ( s<sMin )
    { 
        cout << "Chosen value for the value of s=" << s 
             << " (the center of mass energy) is smaller then the energy threshold, smin=" << sMin
             << ", for the scattering process " << toString(process) << "! Returning 0!\n"; 
        return 0.0; 
    }

    //calculate outgoing particle energies
    double E3 = energyOutParticle3(s, m3, m4);
    double E4 = energyOutParticle4(s, m3, m4);

    //calculate minimum and maximum values of t
    double tMin = tChannelMin(s, m1, m2, m3, m4);
    double tMax = tChannelMax(s, m1, m2, m3, m4);

    int integrationWorkspace = 1000;
    Integration1DimGSLQAGS crossSectionIntegral(tMin/s, tMax/s, &aux, crossSectionProcess12To34Integrand, crossSecIntPrecision, crossSecIntPrecision, integrationWorkspace);

    double crossSection = crossSectionIntegral.evaluate();
    crossSection = s*crossSection;

    //include the Fermi blocking factor for the final states 
    crossSection = crossSection*( 1.0 - fermiDistribution(T, E3 - cP3) )*( 1.0 - fermiDistribution(T, E4 - cP4) );

    //for processes involving completely identical particles, we have to remove extra counting
    if ( process==UUUU || process==DDDD || process==SSSS ||
         process==UBarUBarUBarUBar || process==DBarDBarDBarDBar || process==SBarSBarSBarSBar )
    { 
        crossSection = 0.5*crossSection; 
    }

    return crossSection;
}


void evaluateCrossSectionProcess12To34ToFile(SU3NJL3DCutoffParameters parametersNJL, double T, 
                                             double effChemPotU, double effChemPotD, double effChemPotS, 
                                             double effMassU, double effMassD, double effMassS, 
                                             double propIntPrecision, scatteringProcess process, 
                                             bool largeAngleScatteringContribution, double crossSecIntPrecision,
                                             int numberOfPoints)
{   
    //create arrays to save the calculation of the cross sections as a function of the center of mass energy
    vector<double> sqrtCenterOfMassEnergy(numberOfPoints);
    vector<double> crossSection(numberOfPoints);


    //set incoming and outgoing masses
    double m1, m2, m3, m4;
    inOutMassesGivenScatteringProcess(effMassU, effMassD, effMassS, process, m1, m2, m3, m4);


    //calculate quantities for iteration
    double sqrtCenterOfMassMin = sqrt( centerOfMassEnergyThreshold( m1, m2, m3, m4) + 1E-8 );
    double sqrtCenterOfMassMax = sqrt( sMaximumKlevansky(parametersNJL.getThreeMomentumCutoff(), effMassU, effMassD, effMassS) );
    double sqrtDelta = (sqrtCenterOfMassMax - sqrtCenterOfMassMin)/(numberOfPoints-1);


    int numberThreads = omp_get_max_threads() - 2; //leave one thread free
    cout << "Number of threads being used: " << numberThreads << "\n";
    #pragma omp parallel for schedule(dynamic) num_threads( numberThreads )
    for (int i = 0; i < numberOfPoints; ++i)
    {   
        sqrtCenterOfMassEnergy[i] = sqrtCenterOfMassMin + i*sqrtDelta;
        double s = pow(sqrtCenterOfMassEnergy[i],2);
        crossSection[i] = crossSectionProcess12To34(parametersNJL, T, 
                                                    effChemPotU, effChemPotD, effChemPotS, 
                                                    effMassU, effMassD, effMassS, 
                                                    s, propIntPrecision, process, 
                                                    largeAngleScatteringContribution, crossSecIntPrecision);
        crossSection[i] = inverseGeVSquaredToMiliBarn*crossSection[i];
    }


    //create filename
    string fileName = "crossSection" 
                    + toString(process) 
                    + "_T" + to_string(T) 
                    + "_CPU" + to_string(effChemPotU)
                    + "_CPD" + to_string(effChemPotD)
                    + "_CPS" + to_string(effChemPotS);
    //std::replace( fileName.begin(), fileName.end(), '.', ','); 
    fileName = fileName + ".dat";


    std::ofstream fileCrossSection;
    fileCrossSection.open(fileName, std::ofstream::out | std::ios::trunc);
    fileCrossSection.precision(15);

    fileCrossSection.precision(15);
    fileCrossSection.width(25);   fileCrossSection << "sqrt(s) [GeV]"; 
    fileCrossSection.width(25);   fileCrossSection << "sigma [mb]"; 
    fileCrossSection << endl;

    for (int i = 0; i <numberOfPoints; ++i)
    {
        fileCrossSection.width(25);   fileCrossSection << sqrtCenterOfMassEnergy[i];
        fileCrossSection.width(25);   fileCrossSection << crossSection[i];
        fileCrossSection << std::endl;
    }
    fileCrossSection.close();

}


void evaluateCrossSectionsKlevanskyPaper(SU3NJL3DCutoffParameters parametersNJL, double T, 
                                         double effChemPotU, double effChemPotD, double effChemPotS, 
                                         double effMassU, double effMassD, double effMassS, 
                                         double propIntPrecision, 
                                         bool largeAngleScatteringContribution, double crossSecIntPrecision,
                                         int numberOfPoints)
{   
    ////////////////////////////////////////////////////////////
    // Mu=Md, Cpu_u=Cp_d=Cp_s=0.0
    // uu->uu, ud->ud, us->us, ss->ss
    // uubar->uubar, uubar->ddbar, uubar->ssbar, udbar->udbar, usbar->usbar, ssbar->uubar, ssbar->ssbar

    vector<scatteringProcess> processes = { UDUD, USUS, UUUU, SSSS,
                                            UDBarUDBar, USBarUSBar,
                                            UUBarUUBar, UUBarDDBar, UUBarSSBar,
                                            SSBarUUBar, SSBarSSBar };

    for (int i = 0; i < int(processes.size()); ++i)
    {
        cout << "Calculating cross section for the process: " << toString(processes[i]) << "\n";
        evaluateCrossSectionProcess12To34ToFile(parametersNJL, T, 
                                                effChemPotU, effChemPotD, effChemPotS, 
                                                effMassU, effMassD, effMassS, 
                                                propIntPrecision, processes[i], 
                                                largeAngleScatteringContribution, crossSecIntPrecision,
                                                numberOfPoints);   
    }

}


void evaluateCrossSectionsEqualLightMassesEqualChemicalPotential(
    SU3NJL3DCutoffParameters parametersNJL, double T,                                                    
    double effChemPotU, double effChemPotD, double effChemPotS, 
    double effMassU, double effMassD, double effMassS, 
    double propIntPrecision, 
    bool largeAngleScatteringContribution, double crossSecIntPrecision,
    int numberOfPoints
)
{   
    ////////////////////////////////////////////////////////////
    // Mu=Md, Cpu_u=Cp_d=Cp_s
    // uu->uu, ud->ud, us->us, ss->ss, 
    // uubar->uubar, uubar->ddbar, uubar->ssbar, udbar->udbar ,usbar->usbar, subar->subar, ssbar->uubar, ssbar->ssbar
    // ubarubar->ubarubar, sbarsbar->sbarsbar, ubardbar->ubardbar, ubarsbar->ubarsbar

    vector<scatteringProcess> processes = { UDUD, USUS, UUUU, SSSS,
                                            UDBarUDBar, USBarUSBar,
                                            UUBarUUBar, UUBarDDBar, UUBarSSBar,
                                            SSBarUUBar, SSBarSSBar,
                                            SUBarSUBar, UBarUBarUBarUBar, SBarSBarSBarSBar,
                                            UBarDBarUBarDBar, UBarSBarUBarSBar };

    for (int i = 0; i < int(processes.size()); ++i)
    {
        cout << "Calculating cross section for the process: " << toString(processes[i]) << "\n";
        evaluateCrossSectionProcess12To34ToFile(parametersNJL, T, 
                                                effChemPotU, effChemPotD, effChemPotS, 
                                                effMassU, effMassD, effMassS, 
                                                propIntPrecision, processes[i], 
                                                largeAngleScatteringContribution, crossSecIntPrecision,
                                                numberOfPoints);   
    }

}

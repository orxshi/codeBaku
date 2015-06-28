#ifndef COEFFS_H
#define	COEFFS_H

#include "../Grid/Grid.h"

struct Coeffs
{
    double rhoRef;
    double MachRef;
    double MachAirfoil;
    double MachAir;
    double totalArea;
    double dynPres;
    double liftCoef;
    
    bool dynPresSet;
    bool totalAreaSet;
    
    ofstream out;
    
    Coeffs (const Grid& gr, double rhoRef, double MachAir, double MachAirfoil);
    void getCoeffs (const Grid& gr);
};

#endif	/* COEFFS_H */


#include "Coeffs.h"

Coeffs::Coeffs (const Grid& gr, double rhoRef, double MachAir, double MachAirfoil)
{
    this->rhoRef = rhoRef;    
    this->MachAirfoil = MachAirfoil;
    this->MachAir = MachAir;
    totalArea = 0.;
    dynPres = 0.;    
    dynPresSet = false;
    totalAreaSet = false;
    
    string dir = gr.outputDir;
    string temps = "lift.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);
    out.open (dir);
}

void Coeffs::getCoeffs (const Grid& gr)
{
    if (!totalAreaSet)
    {
        for (const Cell& cll: gr.cell)
        {
            if (cll.bc == BC::SLIP_WALL)
            {
                totalArea += mag (gr.face[cll.face[0]].area);
            }
        }
        
        totalAreaSet = true;
    }
    
    if (!dynPresSet)
    {
        // includes area too!
        
        MachRef = MachAir - MachAirfoil;
        dynPres = 0.5 * rhoRef * pow(MachRef,2.) * totalArea;
        
        dynPresSet = true;
    }
    
    CVector F;
    F[0] = 0.;
    F[1] = 0.;
    F[2] = 0.;
    
    // total force on airfoil
    for (const Cell& cll: gr.cell)
    {
        if (cll.bc == BC::SLIP_WALL)
        {
            F += cll.prim[4] * gr.face[cll.face[0]].area;
        }
    }
    
    liftCoef = F[1] / dynPres;
}




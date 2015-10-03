#include <petscsys.h>
#include "../../../../src/Grid/Grid.h"
#include "../../../../src/Output/Output.h"
#include "../../../../src/Solid/Solid.h"
#include "../../../../src/Coeffs/Coeffs.h"
#include "../../../../src/Solver/Solver.h"
#include "../../../../src/AFT/AFT.h"
#include "Init/Init.h"
#include "MovingGrid/OscAirfoil.h"
#include "MovingGrid/StraightMovingAirfoil.h"
#include "Output/Output.h"

int main(int argc, char** argv)
{
    //PetscInitialize (&argc, &argv, NULL, NULL);
    
    Watch watchSteady;
    Watch watchOscAirfoil;
    
    /*PetscMPIInt rank, n_procs;
    MPI_Comm world = PETSC_COMM_WORLD;
    MPI_Comm_rank (world, &rank);
    MPI_Comm_size (world, &n_procs);*/
    
    string mainDir = createOutputDir();
    
    // background grid
    Grid bg (mainDir, 0);
    bg.read_grid();
    bg.set_grid();
    
    // airfoil grid
    Grid ag (mainDir, 1);
    ag.read_grid();
    ag.set_grid();
    
    OscInit oscInit;
    oscInit.read();    
    oscInit.init (ag);
    oscInit.init (bg);
    
    // push grids to vector
    vector<Grid> grs;    
    grs.push_back(bg);
    grs.push_back(ag);
    
    grs[0].setWallDistance(2);
    grs[1].setWallDistance(3);
    
    //finalGrid.readInput();
    //finalGrid.leastSquaresCoeffs();
    //oscInit.init (finalGrid);

    //Solver solSteady (finalGrid, "SOLVER-STEADY");    
    //solSteady.read ("Solver/solSteady.dat");
    
    ///Solver solOscAirfoil (ag, "SOLVER-OSC-AIRFOIL");
    //solOscAirfoil.read ("Solver/solOscAirfoil.dat");

    // solve steady state
    //SMAirfoil sma (solSteady.dt);
    OscAirfoil oa (1.);
    //sma.read ("MovingGrid/smAirfoil.dat");
    oa.read ("MovingGrid/oscAirfoil.dat");

    //Coeffs coeffs (gr, oscInit.rhoInf, oscInit.Mach, oa.MachAirfoil);
    
    /*sma.getAllFaceVelocities (finalGrid);    
    watchSteady.start();
    (solSteady.implicit) ? solSteady.impl(finalGrid) : solSteady.expl(finalGrid);
    solSteady.petsc.finalize();
    watchSteady.stop();    */
    
    int countr = 0;
    watchOscAirfoil.start();
    
    for (double time=0.; time<0.; time+=1.)
    {
        oa.setAngles (time);
        
        grs[0].cellADT.build (grs[0]);
        grs[1].cellADT.build (grs[1]);        
        grs[0].identifyIBlank (grs[1]);
        grs[1].identifyIBlank (grs[0]);
        
        grs[0].outAllVTK (countr);
        grs[1].outAllVTK (countr);
        
        if (time != 7)
        {
            oa.moveGrid (grs[1]);
        }
    }
    
    
    
    // solve osc airfoil
    for (double time=0.; time<1.; time+=1.)
    {
        cout << "time = " << time << endl;
    
        oa.setAngles (time);
        
        grs[0].cellADT.build (grs[0]);
        grs[1].cellADT.build (grs[1]);        
        grs[0].identifyIBlank (grs[1]);
        grs[1].identifyIBlank (grs[0]);
        
        grs[0].outAllVTK (countr);
        grs[1].outAllVTK (countr);
        
        Grid finalGrid (mainDir, 3);
        
        AFT::aft (grs, finalGrid);
                
        //oa.setAngles (solOscAirfoil.time);
        //oa.getAllFaceVelocities (gr);
        //(solOscAirfoil.implicit) ? solOscAirfoil.impl(gr) : solOscAirfoil.expl(gr);
        //coeffs.getCoeffs (gr);
        //outLiftCoef (coeffs, oa.alpha, solOscAirfoil.time);
        finalGrid.outAllVTK (countr);
        oa.moveGrid (grs[1]);
        
        ++countr;
    }
    watchOscAirfoil.stop();    
    
    /*if (rank == MASTER_RANK)
    {
        //gr.outAllTecplot();
        finalGrid.outAllVTK (0);
        //coeffs.out.close();
        log (mainDir, watchSteady.elapsedTime, "elapsedTimeSteady", watchSteady.unit);
        //log (mainDir, watchOscAirfoil.elapsedTime, "elapsedTimeOscAirfoil", watchOscAirfoil.unit);
        solSteady.log (finalGrid.logDir);
        //solOscAirfoil.log (gr.logDir);
        sma.log (finalGrid.logDir);
        //oa.log (gr.logDir);
    }*/
    
    //PetscFinalize();

    return 0;
}
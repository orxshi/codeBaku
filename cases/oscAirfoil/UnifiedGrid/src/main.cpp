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
    
    Grid bg (mainDir, 0);
    bg.read_grid();
    bg.set_grid();
    
    Grid ag (mainDir, 1);
    ag.read_grid();
    ag.set_grid();
    
    Grid finalGrid (mainDir, 3);
    
    OscInit oscInit;
    oscInit.read();    
    oscInit.init (ag);
    oscInit.init (bg);
    
    vector<Grid> grs;    
    grs.push_back(bg);
    grs.push_back(ag);
    
    Solid solid0;
    Solid solid1;
    
    solid0.phys = 2;
    solid1.phys = 3;
    
    grs[0].setWallDistance(solid0);
    grs[1].setWallDistance(solid1);
        
    for (int i=0; i<grs.size(); ++i)
    {
        cout << ">>> Building grid[" << i << "].cellADT" << endl;
        
        grs[i].cellADT.build (grs[i]);

        cout << ">>> Built grid[" << i << "].cellADT" << endl;
    }

    cout << ">>> Identifying (non)active cells" << endl;

    for (int i=0; i<grs.size(); ++i)
    {
        for (int j=0; j<grs.size(); ++j)
        {
            if (j != i)
            {
                grs[i].identifyIBlank (grs[j]);
            }
        }
    }
    
    
    grs[0].outAllVTK(0);
    grs[1].outAllVTK(0);
    
    
    
    
    AFT::aft (grs, finalGrid);    
    exit(-2);
    
    /*Grid gr (mainDir, 0);
    OscInit oscInit;
    
    gr.read_grid();
    gr.set_grid();

    Solver solSteady (gr, "SOLVER-STEADY");
    solSteady.read ("Solver/solSteady.dat");
    
    Solver solOscAirfoil (gr, "SOLVER-OSC-AIRFOIL");
    solOscAirfoil.read ("Solver/solOscAirfoil.dat");

    oscInit.read();
    oscInit.init (gr);

    // solve steady state
    SMAirfoil sma (solSteady.dt);
    OscAirfoil oa (solOscAirfoil.dt);
    sma.read ("MovingGrid/smAirfoil.dat");
    oa.read ("MovingGrid/oscAirfoil.dat");

    Coeffs coeffs (gr, oscInit.rhoInf, oscInit.Mach, oa.MachAirfoil);
    
    sma.getAllFaceVelocities (gr);
    watchSteady.start();
    (solSteady.implicit) ? solSteady.impl(gr) : solSteady.expl(gr);
    solSteady.petsc.finalize();
    watchSteady.stop();    
    
    int countr = 0;
    watchOscAirfoil.start();
    
    // solve osc airfoil
    for (solOscAirfoil.time=0.; solOscAirfoil.time<solOscAirfoil.finalTime; solOscAirfoil.time+=solOscAirfoil.dt)
    {
        oa.setAngles (solOscAirfoil.time);
        oa.getAllFaceVelocities (gr);
        (solOscAirfoil.implicit) ? solOscAirfoil.impl(gr) : solOscAirfoil.expl(gr);
        coeffs.getCoeffs (gr);
        outLiftCoef (coeffs, oa.alpha, solOscAirfoil.time);
        gr.outAllVTK (countr);
        oa.moveGrid (gr);
        
        ++countr;
    }
    watchOscAirfoil.stop();    
    
    if (rank == MASTER_RANK)
    {
        gr.outAllTecplot();
        gr.outAllVTK (0);
        coeffs.out.close();
        log (mainDir, watchSteady.elapsedTime, "elapsedTimeSteady", watchSteady.unit);
        log (mainDir, watchOscAirfoil.elapsedTime, "elapsedTimeOscAirfoil", watchOscAirfoil.unit);
        solSteady.log (gr.logDir);
        solOscAirfoil.log (gr.logDir);
        sma.log (gr.logDir);
        oa.log (gr.logDir);
    }
    
    PetscFinalize();*/

    return 0;
}




















/*#include "Method/Method.h"

void customUnsteady (Grid& gr)
{
	presetSingle();

	gr.getMMVariables();
    gr.getMovingFaceVelocity();
    gr.apply_BCs();
    
    gr.alpha = 0.016;
    gr.delAlpha = 0.;
    
    SSG (gr);

    gr.steady = false;
    gr.useCFL = false;
    gr.kc = 0.0814;
    gr.alphaMax = 2.51;
    gr.alphaMean = 0.016;    
    
    double dt = gr.dt;
    
    for (double time=0.; time<52; time += dt)
    {
       	gr.time = time;        

       	gr.getMMVariables();
       	gr.getMovingFaceVelocity();
       	gr.apply_BCs();        
        
		SSG (gr);        

        gr.coeffs();
        gr.rotateMesh();
        
        gr.outAllUnsteady();
    }
}

int main(int argc, char** argv)
{
	Grid gr;

    //cout << string(SPLITTER_LEN, '=') << endl;
	
	double startTime = getWallTime();
    customUnsteady (gr);
    double endTime = getWallTime();
    getElapsedTime (startTime, endTime, gr.elapsedTime, gr.elapsedTimeUnit);
    postset (gr);

    return 0;
}*/

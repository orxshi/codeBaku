#include "Solver.h"

void Solver::impl (Grid& gr)
{   
    //Watch wt;
    
    //preSolverCheck (gr);
    
    string dir = gr.outputDir;
    string temps = "res.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);
    
    for (nTimeStep=0; nTimeStep<maxTimeStep; ++nTimeStep)
    {        
        interflux(gr);
        
        switch (linearSolverType)
        {
            case 1:
                gauss_seidel (gr);
                break;
            case 2:
                //MPI_Barrier (PETSC_COMM_WORLD);
                //wt.start();                
                petsc.solveAxb (gr);                
                //FILE *fp;
                //fp=fopen("../out/petscLog", "w");
                //PetscMallocDump(fp);
                //PetscLogDouble pld;
                //PetscMallocGetCurrentUsage (&pld);
                //cout << "pld = " << pld << endl;
                //wt.stop();
                //PetscPrintf (this->petsc.world, "%f\n", wt.elapsedTime);
                //MPI_Barrier (PETSC_COMM_WORLD);
                break;
            default:
                cout << "undefined linear solver in Solver::impl(...)" << endl;
                exit(-2);
                break;
        }
        
        diff_to_cons_prim (gr);
        set_residual (gr);
        
        if (petsc.rank == MASTER_RANK) { outRes(gr.outputDir); }
        gr.apply_BCs();
        
        if (verbose && petsc.rank == MASTER_RANK)
        {
            cout << left << setw(10) << fixed << time;
            cout << setw(10) << nTimeStep;
            cout << scientific << aveRes << endl;
        }        
        
        if (fabs(aveRes) < tol) { break; }
    }
    
    for (Cell& cll: gr.cell)
    {
        cll.oldold_cons = cll.old_cons;
        cll.old_cons = cll.cons;
    }
    
    glo_nTimeStep += nTimeStep;
    ++nImplicitCalls;
}

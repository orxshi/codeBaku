/* 
 * File:   Solver.h
 * Author: Orhan Shibliyev
 *
 * Created on October 2, 2014, 10:00 PM
 */

#ifndef SOLVER_H
#define	SOLVER_H

//#include <petscsys.h>
#include <petscksp.h>
#include "../Grid/Grid.h"
#include "../Output/Output.h"
#include "../InterFlux/Roe/Roe.h"
#include "../Time/Time.h"

//enum class tsOrder_t {UNDEFINED=0, FIRST=1, SECOND=2}; // temporal and spatial orders
//enum class LinearSol_t {UNDEFINED=0, MYGS=1, PETSC=2}; // linear solver for implicit method

struct Solver
{
    int nGaussIter;
    int maxTimeStep;
    int nTimeStep;
    int glo_nTimeStep;
    int nImplicitCalls;
    int tOrder;
    int sOrder;
    int linearSolverType;
    double cfl;
    double dt;
    double finalTime;
    double tol;
    double oversetTol;
    double time;
    double aveRes;    
    bool steady;
    bool implicit;
    bool verbose;
    //tsOrder_t tOrder;
    //tsOrder_t sOrder;
    //LinearSol_t linearSolverType;        
    string instanceName;
    
    /*Vec x, b;
    Mat A;
    PetscScalar* dx;
    PetscInt n;
    PetscInt bs;
    
    MPI_Comm world;*/
    
    Solver (Grid& gr, string instanceName);
    
    void expl (Grid& gr);
    void impl (Grid& gr);
    void interflux (Grid& gr);
    void gauss_seidel (Grid& g);
    void updateVars (Grid& gr);
    double setExpRes (Grid& gr);
    //void preSolver(Grid& gr);
    void log (string fileName);
    void outRes (string fileName);
    void preSolverCheck(const Grid& gr);    
    void read (string fileName);
    void petsc (Grid& gr);
    void set_residual (Grid& g);
    void diff_to_cons_prim(Grid& g);
    
    bool cm (string s, ifstream& in);
    template <class T> void cmh (string s, string membs, T& memb, ifstream& in, bool& found) // helper fnc for cm
    {
        if (!found)
        {
            int isEqual = 0;

            if ( s.compare ( membs ) == isEqual )
            {
                in >> memb;
                found = true;
                
                /*cout << s << endl;
                cout << membs << endl;  
                cout << memb << endl;*/
            }
        }
    }
};

#endif	/* SOLVER_H */


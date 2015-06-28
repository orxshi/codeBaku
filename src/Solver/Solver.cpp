#include "Solver.h"

Solver::Solver (Grid& gr, string instanceName)
{
    //default    
    tOrder = 2;
    sOrder = 1;
    linearSolverType = 1; // MYGS
    nGaussIter = 5;
    maxTimeStep = 10000;
    cfl = 5;
    dt = 1.;
    finalTime = 10.;
    tol = 1e-12;
    oversetTol = 1e-15;
    steady = false;
    implicit = true;
    verbose = true;
    
    // initialize
    time = 0.;
    glo_nTimeStep = 0;
    nImplicitCalls = 0;    
    this->instanceName = instanceName;
    
    /*world = PETSC_COMM_WORLD;
    n = gr.n_in_elm;
    bs = N_VAR;
    
    // set x
    VecCreate (world, &x);
    VecSetType (x, VECSEQ);
    VecSetSizes (x, PETSC_DECIDE, n*bs);
    
    // set b
    VecDuplicate (x, &b);
    
    // set A
    MatCreate (world, &A);
    MatSetType (A, MATSEQBAIJ);
    MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE, n*bs, n*bs);
    MatSeqBAIJSetPreallocation (A, bs, 4, NULL);
    //MatCreateBAIJ (world, bs, PETSC_DECIDE, PETSC_DECIDE, n*bs, n*bs, 1, NULL, 3, NULL, &A);
    // specific to pentagonal mesh . change later
    
    dx = (PetscScalar*)malloc (n * bs * sizeof(PetscScalar));*/
}

void Solver::preSolverCheck (const Grid& gr)
{
    string reason = "nothing";
    
    for (const Cell& cll: gr.cell)
    {
        if (cll.iBlank == iBlank_t::UNDEFINED)
        {
            reason = "iBlank is undefined in Solver::preSolverCheck(...)";
            break;
        }
    }
}

bool Solver::cm (string s, ifstream& in)
{    
    bool found = false;
    
    cmh (s, STRINGTIFY(nGaussIter), nGaussIter, in, found);
    cmh (s, STRINGTIFY(maxTimeStep), maxTimeStep, in, found);
    cmh (s, STRINGTIFY(cfl), cfl, in, found);
    cmh (s, STRINGTIFY(dt), dt, in, found);
    cmh (s, STRINGTIFY(finalTime), finalTime, in, found);
    cmh (s, STRINGTIFY(tol), tol, in, found);
    cmh (s, STRINGTIFY(oversetTol), oversetTol, in, found);
    cmh (s, STRINGTIFY(steady), steady, in, found);
    cmh (s, STRINGTIFY(implicit), implicit, in, found);
    cmh (s, STRINGTIFY(verbose), verbose, in, found);
    cmh (s, STRINGTIFY(tOrder), tOrder, in, found);
    cmh (s, STRINGTIFY(sOrder), sOrder, in, found);
    cmh (s, STRINGTIFY(linearSolverType), linearSolverType, in, found);
    
    return found;
}

void Solver::read (string fileName)
{
    string tmps;
    ifstream in;
    in.open (fileName);
    
    if (in.is_open())
    {
        in >> tmps;
        while ( !in.eof() )
        {
            if ( cm (tmps, in) == false )
            {
                cout << "undefined input in Solver::read(...)" << endl;
                exit(-2);
            }

            in >> tmps;
        }
    }
    else
    {
        cout << "could not open file in Solver::read(...)" << endl;
        exit(-2);
    }

    /*in >> tmps; in >> nGaussIter;
    in >> tmps; in >> maxTimeStep;    
    in >> tmps; in >> cfl;
    in >> tmps; in >> dt;
    in >> tmps; in >> finalTime;
    in >> tmps; in >> tol;
    in >> tmps; in >> steady;
    in >> tmps; in >> implicit;    

    in >> tmps; in >> tmpi;
    if (tmpi == 1)
    {
            tOrder = tsOrder_t::FIRST;
    }
    else if (tmpi == 2)
    {
            tOrder = tsOrder_t::SECOND;
    }
    else
    {
            cout << "tOrder is neither 1 nor 2 in Solver::read(...)" << endl;
            exit(-2);
    }

    in >> tmps; in >> tmpi;
    if (tmpi == 1)
    {
            sOrder = tsOrder_t::FIRST;
    }
    else if (tmpi == 2)
    {
            sOrder = tsOrder_t::SECOND;
    }
    else
    {
            cout << "sOrder is neither 1 nor 2 in Solver::read(...)" << endl;
            exit(-2);
    }

    in >> tmps; in >> verbose;
    
    in >> tmps; in >> tmpi;
    if (tmpi == 1)
    {
            linearSolverType = LinearSol_t::MYGS;
    }
    else if (tmpi == 2)
    {
            linearSolverType = LinearSol_t::PETSC;
    }
    else
    {
            cout << "linearSolverType is neither 1 nor 2 in Solver::read(...)" << endl;
            exit(-2);
    }*/

    in.close();
}
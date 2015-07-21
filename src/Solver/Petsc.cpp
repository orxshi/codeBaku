#include "Solver.h"

void Solver::Petsc::solveAxb (Grid& gr)
{
    VecCreate (world, &x);
    VecSetType (x, VECSTANDARD);
    VecSetSizes (x, PETSC_DECIDE, xGlobalSize);
    VecGetLocalSize (x, &xLocalSize);
    VecGetOwnershipRange (x, &vecFirst, &vecLast);

    // set b
    VecDuplicate (x, &b);
    
    // set A
    MatCreate (world, &A);    
    //MatSetType (A, MATSEQBAIJ);
    MatSetType (A, MATMPIAIJ);
    MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE, xGlobalSize, xGlobalSize);
    //MatSeqBAIJSetPreallocation (A, bs, 4, NULL);
    //MatSeqAIJSetPreallocation (A, 4*bs, NULL);    
    MatMPIAIJSetPreallocation (A, 4*bs, NULL, 4*bs, NULL);
    //MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    //MatCreateBAIJ (world, bs, PETSC_DECIDE, PETSC_DECIDE, n*bs, n*bs, 1, NULL, 3, NULL, &A);
    // specific to pentagonal mesh . change later
    MatGetOwnershipRange (A, &first, &last);
    
    double *dx = NULL;
    int nProcs;
    PetscMPIInt rank;
    MPI_Comm_rank (world, &rank);
    MPI_Comm_size (world, &nProcs);
    
    // set values of b
    for (PetscInt gp=first; gp<last; ++gp)
    {
        PetscInt brow = static_cast <int> (floor(gp/bs));
        PetscInt c = brow + gr.n_bou_elm;
        PetscInt i = gp % bs;
        
        Cell& cll = gr.cell[c];
        
        VecSetValue (b, gp, cll.R[i], INSERT_VALUES);
        // u can use vecsetvalues instead
    }
    
    VecAssemblyBegin (b);
    VecAssemblyEnd (b);
    
    // set values of A
    for (PetscInt gp=first; gp<last; ++gp)
    {
        PetscInt brow = static_cast <int> (floor(gp/bs));
        PetscInt gpDiagStart = brow*bs;
        PetscInt c = brow + gr.n_bou_elm;
        PetscInt i = gp % bs;
        
        Cell& cll = gr.cell[c];
        PetscInt idxn[bs];
        double v[bs];
        
        for (int j=0; j<bs; ++j) { idxn[j] = gpDiagStart + j; }        
        for (int j=0; j<bs; ++j) { v[j] = cll.D[i][j]; }
        
        MatSetValues (A, 1, &gp, bs, idxn, v, INSERT_VALUES);
        
        for (int nn=0; nn<cll.nei.size(); ++nn)
        {
            if (cll.nei[nn] >= gr.n_bou_elm)
            {
                PetscInt neiBcol = cll.nei[nn] - gr.n_bou_elm;
                PetscInt gpNeiStart = neiBcol * bs;
                
                for (int j=0; j<bs; ++j) { idxn[j] = gpNeiStart + j; }
                
                Face& f = gr.face[cll.face[nn]];
                if (c == f.nei[0])
                {
                    for (int j=0; j<bs; ++j) { v[j] = f.M[1][i][j]; }
                }
                else
                {
                    for (int j=0; j<bs; ++j) { v[j] = -f.M[0][i][j]; }
                }
                
                MatSetValues (A, 1, &gp, bs, idxn, v, INSERT_VALUES);
            }
        }
    }
    
    MatAssemblyBegin (A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY);
    
    //------------------------------------------------------
    
    KSPCreate (world, &ksp);
    //KSPSetType (ksp, KSPPREONLY);
    KSPGetPC (ksp, &pc);
    //PCSetType (pc, PCLU);
    KSPSetOperators (ksp, A, A);    
    PCSetType (pc, PCSOR);    
    //PCSORSetIterations (pc, 5, 1);
    //PCSORSetSymmetric (pc, SOR_FORWARD_SWEEP);
    KSPSetType (ksp, KSPGMRES);    
    //KSPSetType (ksp, KSPRICHARDSON);
    //KSPSetTolerances (ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
    KSPSetFromOptions (ksp);    
    KSPSolve (ksp, b, x);
    
    //------------------------------------------------------    
    
    int localSizes[nProcs];
    int recvcounts[nProcs];
    int displs[nProcs];    
    displs[0] = 0;
    
    for (int i=0; i<nProcs; ++i)
    {
        recvcounts[i] = 1;
    }
    for (int i=1; i<nProcs; ++i)
    {
        displs[i] = displs[i-1] + recvcounts[i-1];
    }
    
    //MPI_Gatherv (&xLocalSize, 1, MPI_INT, localSizes, recvcounts, displs, MPI_INT, MASTER_RANK, world);
    MPI_Allgatherv (&xLocalSize, 1, MPI_INT, localSizes, recvcounts, displs, MPI_INT, world);
    
    for (int i=0; i<nProcs; ++i)
    {
        recvcounts[i] = localSizes[i];
    }
    for (int i=1; i<nProcs; ++i)
    {
        displs[i] = displs[i-1] + recvcounts[i-1];
    }
    
    VecGetArray (x, &dx);
    //for (PetscInt gp=first; gp<last; ++gp)
    //{
        //DX[gp] = dx[gp-first];
        
        //MPI_Gatherv (dx, localSizes[rank], MPI_DOUBLE, DX, localSizes, displs, MPI_DOUBLE, MASTER_RANK, world);
        MPI_Allgatherv (dx, localSizes[rank], MPI_DOUBLE, DX, localSizes, displs, MPI_DOUBLE, world);
    //}
    VecRestoreArray (x, &dx);
    
    //MPI_Bcast (DX, xGlobalSize, MPI_DOUBLE, MASTER_RANK, world);
    
    //MPI_Barrier (PETSC_COMM_WORLD);    
    //MPI_Bcast (DX + first, xLocalSize, MPI_DOUBLE, rank, world);
    //MPI_Bcast (DX + 0, xLocalSize, MPI_DOUBLE, 0, world);
    /*MPI_Status status;
    if (rank == 0)
    {
        MPI_Send (DX, xLocalSize, MPI_DOUBLE, 1, 0, world);
        MPI_Recv (DX + xLocalSize, xLocalSize, MPI_DOUBLE, 1, 0, world, &status);
    }
    else if (rank == 1)
    {
        MPI_Send (DX + xLocalSize, xLocalSize, MPI_DOUBLE, 0, 0, world);
        MPI_Recv (DX, xLocalSize, MPI_DOUBLE, 0, 0, world, &status);
    }*/
    //MPI_Barrier (PETSC_COMM_WORLD);
    
    for (PetscInt gp=0; gp<xGlobalSize; ++gp)
    {
        PetscInt brow = static_cast <int> (floor(gp/bs));
        PetscInt c = brow + gr.n_bou_elm;
        PetscInt i = gp % bs;        
        Cell& cll = gr.cell[c];
        
        cll.dQ[i] = DX[gp];
    }
    
    finalize();
}

void Solver::Petsc::finalize()
{
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    //KSPDestroy(&ksp);
    //PetscFree (DX);
    //DX = NULL;
}

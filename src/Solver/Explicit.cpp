#include "Solver.h"

void Solver::updateVars (Grid& gr)
{
    for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
    {
        Cell& cll = gr.cell[c];
        
        /*if (useCFL)
        {
            dt = cfl * cll.vol / cll.sigma;
        }*/
        
        if (tOrder == 1 || nTimeStep < 2)
        {
            for (int i=0; i<N_VAR; ++i)
            {
                cll.cons[i] += dt * cll.R[i] / cll.vol;
            }
        }
        else if (tOrder == 2)
        {
            for (int i=0; i<N_VAR; ++i)
            {
                cll.cons[i] = ( 2. * dt * cll.R[i] / cll.vol + 4.*cll.old_cons[i] - cll.oldold_cons[i] ) / 3.;
            }
        }
        else
        {
            cout << "time derivative is not 1 or 2 (explicit scheme)" << endl;
            exit(-2);
        }        
        
        cll.cons_to_prim();
    }
}

double Solver::setExpRes (Grid& gr)
{
    Vector<N_VAR> res;
    res.fill(0.);    
    
    // calculate cll.R with updated values. updateVars should be called before this.
    for (Face& f: gr.face)
    {
        roeflx(f, gr.cell[f.nei[0]], gr.cell[f.nei[1]]);
    }
    
    if (tOrder == 1 || nTimeStep < 2)
    {
        if (steady)
        {
            for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
            {
                Cell& cll = gr.cell[c];

                for (int i=0; i<N_VAR; ++i)
                {
                    res[i] += pow( cll.cons[i] - cll.old_cons[i], 2 );
                }
            }
        }
        else
        {
            for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
            {
                Cell& cll = gr.cell[c];

                for (int i=0; i<N_VAR; ++i)
                {
                    double LHS = (cll.cons[i] - cll.old_cons[i]) * cll.vol / dt;
                    double RHS = cll.R[i];

                    res[i] += pow(LHS - RHS, 2);
                }
            }
        }
    }
    else if (tOrder == 2)
    {
        if (steady)
        {
            for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
            {
                Cell& cll = gr.cell[c];

                for (int i=0; i<N_VAR; ++i)
                {
                    res[i] += pow( 3.*cll.cons[i] - 4.*cll.old_cons[i] + cll.oldold_cons[i], 2 );
                }
            }
        }
        else
        {
            for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
            {
                Cell& cll = gr.cell[c];

                for (int i=0; i<N_VAR; ++i)
                {
                    double LHS = 0.5 * (3.*cll.cons[i] - 4.*cll.old_cons[i] + cll.oldold_cons[i]) * cll.vol / dt;
                    double RHS = cll.R[i];

                    res[i] += pow(LHS - RHS, 2);
                }
            }
        }
    }
    else
    {
        cout << "time derivative is not 1 or 2 in Solver::setExpRes(...)" << endl;
        exit(-2);
    }
    

    //for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
    //{
    //    Cell& cll = gr.cell[c];
    //    
    //    for (int i=0; i<N_VAR; ++i)
    //    {
    //        resOuter[i] += pow( cll.cons[i] - cll.old_cons[i], 2 );
    //    }
    ///

    for (int i=0; i<N_VAR; ++i)
    {
        //resOuter[i] /= cell.size();
        res[i] /= (gr.n_in_elm);
        res[i] = sqrt(res[i]);
    }

    double averageRes = 0.;
    for (int i=0; i<N_VAR; ++i)
    {
        averageRes += res[i];
    }
    averageRes /= N_VAR;
    
    if ( isnan(averageRes) ) { cout << "averageRes is NAN in setExpRes()" << endl; exit(-2); }
    if ( isinf(averageRes) ) { cout << "averageRes is INF in setExpRes()" << endl; exit(-2); }

    return averageRes;
}

void Solver::expl (Grid& gr)
{
    preSolverCheck (gr);
    
    string dir = gr.outputDir;
    string temps = "res.dat";
    string slash = "/";
    dir.append (slash);
    dir.append (temps);
    
    double aveRes;
    
    for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
    {
        gr.cell[c].oldold_cons.fill(0.);
        gr.cell[c].old_cons.fill(0.);
    }
    
    for (nTimeStep=0; nTimeStep<maxTimeStep; ++nTimeStep)
    {
        for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
        {
            Cell& cll = gr.cell[c];
            
            cll.sigma = 0.;
            cll.R.fill(0.);
        }
        
        if (sOrder == 2) { gr.leastSquaresGrad(); }
        
        for (Face& f: gr.face)
        {
            roeflx(f, gr.cell[f.nei[0]], gr.cell[f.nei[1]]);
        }
         
        updateVars(gr);
        gr.apply_BCs();        
        
        aveRes = setExpRes(gr);
        
        if (verbose)
        {
            cout << left << setw(10) << time;
            cout << setw(10) << nTimeStep;
            //if (!useCFL) { cout << setw(10) << fixed << setprecision(3) << time; }
            cout << scientific << aveRes << endl;
        }
        
        outRes(gr.outputDir);
        
        if (aveRes < tol)
        //if (steady && fabs(aveRes) < tol)
        {
            break;
        }
        
        for (int c=gr.n_bou_elm; c<gr.cell.size(); ++c)
        {
            Cell& cll = gr.cell[c];
            
            cll.oldold_cons = cll.old_cons;
            cll.old_cons = cll.cons;
        }
        
        //if (!useCFL) { time += dt; }
        time += dt;
    }
}

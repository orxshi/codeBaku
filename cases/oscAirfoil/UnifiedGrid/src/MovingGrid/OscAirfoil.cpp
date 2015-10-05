#include "OscAirfoil.h"

void OscAirfoil::setAngles(const double time)
{
    double tmp = 2. * kc * fabs(MachAirfoil);
    //alpha    = alphaMean + alphaMax * sin( tmp * (time+dt) ); // degree
    //delAlpha = -dt * tmp * alphaMax * cos( tmp * (time+dt) ); // degree
    
    alpha    = alphaMean + alphaMax * sin( tmp * time ); // degree
    delAlpha = -tmp * alphaMax * cos( tmp * time ); // degree
}

void OscAirfoil::moveEdge()
{
    CVector r;
    
    double delAlphaRad = delAlpha * DEG_TO_RAD;
    double alphaRad = alpha * DEG_TO_RAD;
    
    r = vertex1 - centerAirfoil;
    rotateVectorAroundPoint2D (centerAirfoil, delAlphaRad, r, vertex2);

    r = vertex0 - centerAirfoil;
    rotateVectorAroundPoint2D (centerAirfoil, delAlphaRad, r, vertex3);

    CVector vel;
    vel[0] = MachAirfoil * cos (alphaRad);
    vel[1] = MachAirfoil * sin (alphaRad);
    vel[2] = 0.;

    displacePoint2D (vertex2, vel, dt, vertex2);
    displacePoint2D (vertex3, vel, dt, vertex3);
}

void OscAirfoil::moveGrid (Grid& gr)
{
    CVector r;

    double dar = delAlpha * DEG_TO_RAD;
    
    for (Point& p: gr.pt)
    {
        r = p.dim - centerAirfoil;
        rotateVectorAroundPoint2D (centerAirfoil, dar, r, p.dim);
    }
    
    for (Cell& c: gr.cell) { c.set_centroid (gr.pt); }

    for (Face& f: gr.face)
    {
        f.set_centroid (gr.pt);
        f.set_area (gr.pt);
    }

    for (Cell& c: gr.cell)
    {
        double oldU = c.prim[1];
        double oldV = c.prim[2];

        c.prim[1] = oldU * cos(dar) - oldV * sin(dar);
        c.prim[2] = oldU * sin(dar) + oldV * cos(dar);
        
        c.prim_to_cons();
    }
}

void OscAirfoil::read (string fileName)
{
    string tmps;
    ifstream in;
    in.open (fileName);
    
    if (in.is_open())
    {
        in >> tmps; in >> alphaMean;
        in >> tmps; in >> alphaMax;
        in >> tmps; in >> kc;
        in >> tmps; in >> MachAirfoil;
        in >> tmps; in >> centerAirfoil[0];
        in >> tmps; in >> centerAirfoil[1];
        in >> tmps; in >> centerAirfoil[2];
    }
    else
    {
        cout << "could not open file in OscAirfoil::read(...)" << endl;
        exit(-2);
    }
    
    in.close();
}

void OscAirfoil::log (string fileName)
{
    ofstream out;
    out.open (fileName, ofstream::app);
    
    if (out.is_open())
    {
        out << endl;
        out << "OSCILLATING AIRFOIL" << endl;
        
        out << "alphaMean = " << alphaMean << endl;
        out << "alphaMax = " << alphaMax << endl;
        out << "kc = " << kc << endl;
        out << "MachAirfoil = " << MachAirfoil << endl;
        out << "centerAirfoil[0] = " << centerAirfoil[0] << endl;
        out << "centerAirfoil[1] = " << centerAirfoil[1] << endl;
        out << "centerAirfoil[2] = " << centerAirfoil[2] << endl;
        out << "dt = " << dt << endl;
    }
    else
    {
        cout << "could not open file in OscAirfoil::log(...)" << endl; // log as well
        exit(-2);
    }
    
    out.close();
}

void OscAirfoil::interFromOldTS (Grid& curGrid, Grid& oldGrid)
{
    for (int c=curGrid.n_bou_elm; c<curGrid.cell.size(); ++c)    
    {
        Cell& cll = curGrid.cell[c];
        
        if (cll.belonging == 0) continue; // airfoil or new grid
        
        CVector r, p;

        double dar = delAlpha * DEG_TO_RAD;
        
        p = cll.cnt;
        r = p - centerAirfoil;
        rotateVectorAroundPoint2D (centerAirfoil, -dar, r, p);
        
        ADT::ADTPoint vec;
        
        for (int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = p[i];
            vec.dim[i*2+1] = p[i];
        }       
        
        int ind = oldGrid.cellADT.search (vec);
        
        if (ind != -1)
        {
            Cell& oc = oldGrid.cell[ind];
            cll.prim = oc.prim;
            cll.cons = oc.cons;
        }
        else
        {
            cout << "could not find corresponding cell in Grid::interFromOldTS(...)" << endl;
            cout << "c = " << c << endl;
            cout << "delAlpha = " << delAlpha << endl;
            cout << "p[0] = " << p[0] << endl;
            cout << "p[1] = " << p[1] << endl;
            cout << "p[2] = " << p[2] << endl;
            cout << "cll.belonging = " << cll.belonging << endl;
            cout << "curGrid.n_bou_elm = " << curGrid.n_bou_elm << endl;
            cout << "curGrid.cell.size() = " << curGrid.cell.size() << endl;
            cout << "oldGrid tree size = " << oldGrid.cellADT.idsInTree.size() << endl;
            cout << "idx = " << oldGrid.cellADT.root->p->idx << endl;
            cout << "idx = " << oldGrid.cellADT.root->p->dim[0] << endl;
            cout << "idx = " << oldGrid.cellADT.root->p->dim[1] << endl;
            cout << "idx = " << oldGrid.cellADT.root->p->dim[2] << endl;
            cout << "idx = " << oldGrid.cellADT.root->p->dim[3] << endl;
            cout << "idx = " << oldGrid.cellADT.root->p->dim[4] << endl;
            cout << "idx = " << oldGrid.cellADT.root->p->dim[5] << endl;
            cout << "cc = " << oldGrid.cellADT.root->c[0] << endl;
            cout << "cc = " << oldGrid.cellADT.root->c[1] << endl;
            cout << "cc = " << oldGrid.cellADT.root->c[2] << endl;
            cout << "cc = " << oldGrid.cellADT.root->c[3] << endl;
            cout << "cc = " << oldGrid.cellADT.root->c[4] << endl;
            cout << "cc = " << oldGrid.cellADT.root->d[5] << endl;
            cout << "dd = " << oldGrid.cellADT.root->d[0] << endl;
            cout << "dd = " << oldGrid.cellADT.root->d[1] << endl;
            cout << "dd = " << oldGrid.cellADT.root->d[2] << endl;
            cout << "dd = " << oldGrid.cellADT.root->d[3] << endl;
            cout << "dd = " << oldGrid.cellADT.root->d[4] << endl;
            cout << "dd = " << oldGrid.cellADT.root->d[5] << endl;
            cout << "tt = " << vec.dim[0] << endl;
            cout << "tt = " << vec.dim[1] << endl;
            cout << "tt = " << vec.dim[2] << endl;
            cout << "tt = " << vec.dim[3] << endl;
            cout << "tt = " << vec.dim[4] << endl;
            cout << "tt = " << vec.dim[5] << endl;
            for (int v: cll.vtx)
            {
                cout << "v = " << v << endl;
            }
            for (int v: cll.vtx)
            {
                cout << "curGrid.pt[v].dim = " << curGrid.pt[v].dim[0] << endl;
                cout << "curGrid.pt[v].dim = " << curGrid.pt[v].dim[1] << endl;
                cout << "curGrid.pt[v].dim = " << curGrid.pt[v].dim[2] << endl;
            }
            
            for (int i=0; i<oldGrid.cellADT.idsInTree.size(); ++i)
            {                
                if (oldGrid.cellADT.idsInTree[i] == c)
                {
                    
                    if (oldGrid.cellADT.addrsInTree[i] != NULL)
                    {
                        if (oldGrid.cellADT.addrsInTree[i]->p != NULL)
                        {
                            ADT::Node* node = oldGrid.cellADT.addrsInTree[i];
                            ADT::ADTPoint* ap = node->p;
                        
                            cout << "baku = " << ap->idx << endl;
                            cout << "c = " << node->c[0] << endl;
                            cout << "c = " << node->c[1] << endl;
                            cout << "c = " << node->c[2] << endl;
                            cout << "c = " << node->c[4] << endl;
                            cout << "c = " << node->c[4] << endl;
                            cout << "c = " << node->c[5] << endl;
                            
                            cout << "d = " << node->d[0] << endl;
                            cout << "d = " << node->d[1] << endl;
                            cout << "d = " << node->d[2] << endl;
                            cout << "d = " << node->d[3] << endl;
                            cout << "d = " << node->d[4] << endl;
                            cout << "d = " << node->d[5] << endl;
                            
                            bool cubesOverlap = oldGrid.cellADT.doCubesOverlap (node, vec);
                            bool cmp = oldGrid.cellADT.compareFunction (node, vec);
                            
                            cout << "cubesOverlap = " << cubesOverlap << endl;
                            cout << "cmp = " << cmp << endl;
                            cout << "cmp = " << node->p->vertices.size() << endl;
                            cout << "cmp = " << vec.dim[0] << endl;
                            cout << "cmp = " << vec.dim[1] << endl;
                            cout << "cmp = " << vec.dim[2] << endl;
                            cout << "cmp = " << vec.dim[3] << endl;
                            cout << "cmp = " << vec.dim[4] << endl;
                            cout << "cmp = " << vec.dim[5] << endl;
                            
                            exit(-2);
                        }
                    }
                }
            }
            
            curGrid.outAllVTK (1);
            /*for (int i: oldGrid.cellADT.idsInTree)
            {
                cout << "ids = " << i << endl;
            }*/
            exit(-2);
        }
    }
}
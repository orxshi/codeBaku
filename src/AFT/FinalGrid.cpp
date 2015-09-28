#include "AFT.h"

namespace AFT
{
    void addToPointADT (const Point& p, PointADT& pointADT, int index)
    {
        ADT::ADTPoint vec = pointADT.createADTPoint (p.dim, p.dim);
        vec.idx = index;
        bool tempBool;
        pointADT.insert (vec, pointADT.root, tempBool);
    }
    
    bool cmpPoints (const ::Point& p0, const ::Point& p1)
    {
        int counter = 0;
        for (int d=0; d<N_DIM; ++d)
        {
            if ( fabs(p0.dim[d] - p1.dim[d]) < 1e-10 )
            {
                ++counter;
            }
        }

        if (counter == N_DIM)
        {
            return true;
        }
        
        return false;
    }
    
    bool pointExistsForCreateCells(const ::Point& refPoint, const vector<::Point>& points, int& index)
    {
        for (int ip=0; ip<points.size(); ++ip)
        {
            if ( cmpPoints (refPoint, points[ip]) )
            {
                index = ip;
                return true;
            }
        }
        
        return false;
    }
    
    bool cellExistsForCreateCells (const Cell& refCell, const vector<Cell>& cell, int& index)
    {
        for (int c=0; c<cell.size(); ++c)
        {
            int counter = 0;
            for (int d=0; d<N_DIM; ++d)
            {
                if ( fabs(refCell.cnt[d] - cell[c].cnt[d]) < 0.0001 )
                {
                    ++counter;
                }
            }

            if (counter == N_DIM)
            {
                index = c;
                return true;
            }
        }
        
        return false;
    }
    
    bool faceExistsForCreateCells (const Face& refFace, const vector<Face>& face, int& index)
    {
        for (int c=0; c<face.size(); ++c)
        {
            int counter = 0;
            for (int d=0; d<N_DIM; ++d)
            {
                if ( fabs(refFace.cnt[d] - face[c].cnt[d]) < 0.0001 )
                {
                    ++counter;
                }
            }
            
            /*if (refFace.cnt[0] < -6.98640 && refFace.cnt[0] > -6.98642)
            {
                if (face[c].cnt[0] < -6.98640 && face[c].cnt[0] > -6.98642)
                {
                    cout << refFace.cnt[0] << endl;
                    cout << refFace.cnt[1] << endl;
                    cout << refFace.cnt[2] << endl;
                    
                    cout << face[c].cnt[0] << endl;
                    cout << face[c].cnt[1] << endl;
                    cout << face[c].cnt[2] << endl;
                    
                    cout << counter << endl;
                    
                    for (int d=0; d<N_DIM; ++d)
                    {
                        if (refFace.cnt[d] == face[c].cnt[d])
                        {
                            cout << "equal_" << d << endl;
                        }
                        else
                        {
                            cout << refFace.cnt[d] << endl;
                            cout << face[c].cnt[d] << endl;
                        }
                    }
                    
                    cin.ignore();
                }
            }*/

            if (counter == N_DIM)
            {
                index = c;
                return true;
            }
        }
        
        return false;
    }
    
    bool strictlyInteriorPass (const Face& f, const iBlank_t& LCIBlank, const iBlank_t& RCIBlank,
            const bool LCTrim, const bool RCTrim)
    {
        if ( f.bouType != face_t::INTERIOR )
        {
            return false;
        }
        else if ( f.nei.size() != 2 )
        {
            return false;
        }
        else if (LCIBlank != iBlank_t::FIELD || RCIBlank != iBlank_t::FIELD)
        {
            return false;
        }
        else if (LCTrim != false || RCTrim != false)
        {
            return false;
        }
        
        return true;
    }
    
    bool ghostPass (const Face& f, const iBlank_t& LCIBlank, const bool LCTrim)
    {
        if (f.bouType != face_t::BOUNDARY)
        {
            return false;
        }
        else if (LCIBlank != iBlank_t::FIELD || LCTrim != false)
        {
            return false;
        }
        
        return true;
    }
    
    bool intergridPass (const Face& f, const iBlank_t& LCIBlank, const iBlank_t& RCIBlank,
            const bool LCTrim, const bool RCTrim)
    {
        if (f.bouType != face_t::INTERIOR)
        {
            return false;
        }
        else if (LCIBlank != iBlank_t::FIELD || RCIBlank != iBlank_t::FIELD)
        {
            return false;
        }
        else if (LCTrim == false && RCTrim == false)
        {
            return false;
        }
        else if (LCTrim == true && RCTrim == true)
        {
            return false;
        }
                
        return true;
    }
    
    void addGhostCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<::Point>& pt, Grid& finalGrid, PointADT& fgp, PointADT& fgcc)
    {
        //int nExists = 0;
        const Cell* GC;
        GC = &cell[ f.nei[1] ];
            
        // create new face
        Face newFace = f;
        newFace.vtx.clear();
        newFace.nei.clear();
        finalGrid.face.push_back(newFace);
        int iRefFace = finalGrid.face.size()-1;
        Face& refFace = finalGrid.face.back();
        //-----------------------------------------------------------------------                    
        int iRefC;
        Cell* refC;
        
        Cell newCell = *GC;
        newCell.nei.clear();
        newCell.face.clear();
        finalGrid.cell.push_back(newCell);
        iRefC = finalGrid.cell.size()-1;
        Point cellCent;
        cellCent.dim = newCell.cnt;
        addToPointADT (cellCent, fgcc, iRefC);

        refC = &finalGrid.cell[iRefC];
        refFace.nei.push_back(iRefC);
        refC->face.push_back(iRefFace);
        //-----------------------------------------------------------------------                    

        // loop through vertices of face
        for (int v=0; v<f.vtx.size(); ++v)
        {
            // create new point
            ::Point p = pt[ f.vtx[v] ];
            Point pp;
            pp.dim = p.dim;

            //int pointIndex1;
            int pointIndex;
            //bool ptExists = pointExistsForCreateCells(p, finalGrid.pt, pointIndex);            
            bool ptExists = pointExists (pp, fgp, pointIndex);
                        
            /*if (ptExists != ptExists1)
            {
                cout << "not equal" << endl;
                cout << "ptExists1 = " << ptExists1 << endl;
                cout << "ptExists = " << ptExists << endl;
                
                cout << "p.dim[0] = " << p.dim[0] << endl;
                cout << "p.dim[1] = " << p.dim[1] << endl;
                cout << "p.dim[2] = " << p.dim[2] << endl;
                
                cout << "pointIndex1 = " << pointIndex1 << endl;
                cout << "pointIndex = " << pointIndex << endl;
                exit(-2);
            }*/
            
            //if (f.vtx[v] == 8)
            /*{
                cout << "exists = " << ptExists << endl;
                cout << "pointIndex = " << pointIndex << endl;
                
                cout << "p.dim[0] = " << p.dim[0] << endl;
                cout << "p.dim[1] = " << p.dim[1] << endl;
                cout << "p.dim[2] = " << p.dim[2] << endl;
                
                //cout << "finalGrid.pt[3].dim[0] = " << finalGrid.pt[3].dim[0] << endl;
                //cout << "finalGrid.pt[3].dim[1] = " << finalGrid.pt[3].dim[1] << endl;
                //cout << "finalGrid.pt[3].dim[2] = " << finalGrid.pt[3].dim[2] << endl;
                cout << "finalGrid.pt.size() = " << finalGrid.pt.size() << endl;
                cout << "fgp.sizeTree = " << fgp.sizeTree << endl;
                //exit(-2);
                cin.ignore();
            }*/

            if (!ptExists)
            {
                pp.belonging = finalGrid.id;
                finalGrid.pt.push_back (p);
                pointIndex = finalGrid.pt.size() - 1;
                addToPointADT (pp, fgp, pointIndex);
            }
            /*else
            {
                cout << "exists" << endl;
                cout << "point = " << f.vtx[v] << endl;
                cout << "pointIndex = " << pointIndex << endl;
                
                cout << "p.dim[0] = " << p.dim[0] << endl;
                cout << "p.dim[1] = " << p.dim[1] << endl;
                cout << "p.dim[2] = " << p.dim[2] << endl;
                
                cout << "finalGrid.pt[3].dim[0] = " << finalGrid.pt[3].dim[0] << endl;
                cout << "finalGrid.pt[3].dim[1] = " << finalGrid.pt[3].dim[1] << endl;
                cout << "finalGrid.pt[3].dim[2] = " << finalGrid.pt[3].dim[2] << endl;
                
                exit(-2);
            }*/

            refFace.vtx.push_back( pointIndex );
            
            for (int cv=0; cv<refC->vtx.size(); ++cv)
            {
                bool match = cmpPoints( p, pt[ GC->vtx[cv] ] );

                if (match)
                {
                    refC->vtx[cv] = pointIndex;
                    break;
                }
            }
        }
        
        refC = NULL;        
        GC = NULL;
        
        //cout << "nExists = " << nExists << endl;
        //cin.ignore();
    }
    
    void addCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<::Point>& pt, Grid& finalGrid, PointADT& fgp, PointADT& fgcc)
    {
        int neiSize = f.nei.size();
        const Cell* LRC[ neiSize ];
        for (int i=0; i<neiSize; ++i)
        {
            LRC[i] = &cell[ f.nei[i] ];
        }
            
        // create new face
        Face newFace = f;
        newFace.vtx.clear();
        newFace.nei.clear();
        finalGrid.face.push_back(newFace);
        int iRefFace = finalGrid.face.size()-1;
        //Face& refFace = finalGrid.face.back();
        //-----------------------------------------------------------------------                    
        int iRefC[neiSize];
        bool CExists[neiSize];

        for (int i=0; i<neiSize; ++i)
        {
            //CExists[i] = cellExistsForCreateCells (*LRC[i], finalGrid.cell, iRefC[i]);            
            Point cellCent;
            cellCent.dim = (*LRC[i]).cnt;            
            CExists[i] = pointExists (cellCent, fgcc, iRefC[i]);

            if (!CExists[i])
            {
                Cell newCell = *LRC[i];
                newCell.nei.clear();
                newCell.face.clear();
                finalGrid.cell.push_back(newCell);
                iRefC[i] = finalGrid.cell.size()-1;                
                Point cellCent;
                cellCent.dim = newCell.cnt;
                addToPointADT (cellCent, fgcc, iRefC[i]);
            }
            
            finalGrid.face.back().nei.push_back(iRefC[i]);
            finalGrid.cell[iRefC[i]].face.push_back(iRefFace);
        }
        //-----------------------------------------------------------------------                    
        if (neiSize == 2)
        {
            finalGrid.cell[iRefC[0]].nei.push_back(iRefC[1]);
            finalGrid.cell[iRefC[1]].nei.push_back(iRefC[0]);
        }

        // loop through vertices of face
        for (int v=0; v<f.vtx.size(); ++v)
        {
            // create new point
            ::Point p = pt[ f.vtx[v] ];
            Point pp;
            pp.dim = p.dim;

            int pointIndex = -1;
            //bool ptExists = pointExistsForCreateCells(p, finalGrid.pt, pointIndex);
            bool ptExists = pointExists (pp, fgp, pointIndex);

            if (!ptExists)
            {
                pp.belonging = finalGrid.id;
                finalGrid.pt.push_back (p);
                pointIndex = finalGrid.pt.size() - 1;
                addToPointADT (pp, fgp, pointIndex);
            }

            finalGrid.face.back().vtx.push_back( pointIndex );
        }
                
        for (int i=0; i<neiSize; ++i)
        {
            LRC[i] = NULL;
        }
    }
    
    void modifyCellVertices (Grid& finalGrid, const Grid& newGrid, const vector<Grid>& gr, PointADT& fgp)
    {
        bool match;
        int pointIndex;
        
        for (int ic=finalGrid.n_bou_elm; ic<finalGrid.cell.size(); ++ic)
        {
            Cell& c = finalGrid.cell[ic];
            
            for (int& v: c.vtx)
            {
                if (c.belonging == 0 || c.belonging == 1)
                {
                    Point pp;
                    pp.dim = gr[c.belonging].pt[v].dim;
                    match = pointExists (pp, fgp, pointIndex);
                }
                else if (c.belonging == 2)
                {
                    Point pp;
                    pp.dim = newGrid.pt[v].dim;
                    match = pointExists (pp, fgp, pointIndex);
                }
                else
                {
                    cout << "!!! Error: c.belonging != 0 || c.belonging != 1 || c.belonging != 2 in modifyCellVertices(...)" << endl;
                    cout << "!!! Error: c.belonging = " << c.belonging << endl;
                    exit(-2);
                }
                
                if (match)
                {
                    v = pointIndex;
                }
                else
                {
                    cout << "!!! Error: no match in modifyCellVertices(...)" << endl;
                    cout << "!!! Error: c.belonging = " << c.belonging << endl;
                    cout << "!!! Error: dim[0] = " << gr[c.belonging].pt[v].dim[0] << endl;
                    cout << "!!! Error: dim[1] = " << gr[c.belonging].pt[v].dim[1] << endl;
                    cout << "!!! Error: dim[2] = " << gr[c.belonging].pt[v].dim[2] << endl;
                    exit(-2);
                }
                
                /*for (int ip=0; ip<finalGrid.pt.size(); ++ip)
                {
                    Point& p = finalGrid.pt[ip];
                    
                    match = false;
                    if (c.belonging == 0 || c.belonging == 1)
                    {
                        match = cmpPoints( p, gr[c.belonging].pt[v]  );                        
                    }
                    else if (c.belonging == 2)
                    {
                        match = cmpPoints( p, newGrid.pt[v]  );
                    }
                    else
                    {
                        cout << "!!! Error: c.belonging != 0 || c.belonging != 1 || c.belonging != 2 in modifyCellVertices(...)" << endl;
                        cout << "!!! Error: c.belonging = " << c.belonging << endl;
                        exit(-2);
                    }

                    if (match)
                    {
                        v = ip;
                        break;
                    }
                }
                
                if (match == false)
                {
                    cout << "!!! Error: no match in modifyCellVertices(...)" << endl;
                    cout << "!!! Error: c.belonging = " << c.belonging << endl;
                    cout << "!!! Error: dim[0] = " << gr[c.belonging].pt[v].dim[0] << endl;
                    cout << "!!! Error: dim[1] = " << gr[c.belonging].pt[v].dim[1] << endl;
                    cout << "!!! Error: dim[2] = " << gr[c.belonging].pt[v].dim[2] << endl;
                    exit(-2);
                }*/
            }
            
            c.belonging = finalGrid.id;
        }
    }
    
    void addIntergridCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<::Point>& pt, Grid& finalGrid, const Grid& newGrid, PointADT& fgp, PointADT& fgcc, PointADT& fgfc)
    {
        int neiSize = f.nei.size();
        const Cell* LRC[ neiSize ];
        for (int i=0; i<neiSize; ++i)
        {
            LRC[i] = &cell[ f.nei[i] ];
        }
            
        // create new face
        Face newFace = f;
        newFace.vtx.clear();
        newFace.nei.clear();
        finalGrid.face.push_back(newFace);
        int iRefFace = finalGrid.face.size()-1;
        Face& refFace = finalGrid.face.back();
        //-----------------------------------------------------------------------                    
        int iRefC[neiSize];
        bool CExists[neiSize];

        for (int i=0; i<neiSize; ++i)
        {
            //CExists[i] = cellExistsForCreateCells (*LRC[i], finalGrid.cell, iRefC[i]);
            Point cellCent;
            cellCent.dim = (*LRC[i]).cnt;            
            CExists[i] = pointExists (cellCent, fgcc, iRefC[i]);

            if (!CExists[i])
            {
                int faceIndex;
                //bool existsInNewgrid = faceExistsForCreateCells (newFace, newGrid.face, faceIndex);
                //bool existsInNewgrid = cellExistsForCreateCells (*LRC[i], newGrid.cell, iRefC[i]);
                Point faceCent;
                faceCent.dim = newFace.cnt;
                bool existsInNewgrid = pointExists (faceCent, fgfc, faceIndex);
                
                Cell newCell;                
                if (!existsInNewgrid)
                {
                    cout << "!existsInNewgrid in addIntergridCell(...)" << endl;
                    exit(-2);
                    //newCell = *LRC[i];
                }
                else
                {
                    //bool tmpBool = cellExistsForCreateCells (newGrid.cell[newGrid.face[faceIndex].nei[0]], finalGrid.cell, iRefC[i]);                    
                    Point cellCent;
                    cellCent.dim = newGrid.cell[newGrid.face[faceIndex].nei[0]].cnt;
                    bool tmpBool = pointExists (cellCent, fgcc, iRefC[i]);
                    
                    if (tmpBool)
                    {
                        goto honolulu;
                    }
                    else
                    {
                        if (newGrid.face[faceIndex].nei.size() == 1)
                        {
                            iRefC[i] = newGrid.face[faceIndex].nei[0];
                        }
                        else
                        {
                            cout << "newGrid.face[faceIndex].nei.size() != 1" << endl;
                            cout << "newGrid.face[faceIndex].nei.size() = " << newGrid.face[faceIndex].nei.size() << endl;
                            exit(-2);
                        }

                        newCell = newGrid.cell[ iRefC[i] ];
                    }
                }
                
                newCell.nei.clear();
                newCell.face.clear();
                finalGrid.cell.push_back(newCell);
                iRefC[i] = finalGrid.cell.size()-1;                
                Point cellCent;
                cellCent.dim = newCell.cnt;
                addToPointADT (cellCent, fgcc, iRefC[i]);
            }
            
            honolulu:;
            refFace.nei.push_back(iRefC[i]);
            finalGrid.cell[iRefC[i]].face.push_back(iRefFace);
        }
        //-----------------------------------------------------------------------
        
        /*if (CExists[0] == true && CExists[1] == true)
        {
            cout << "both are true" << endl;
            cin.ignore();
        }*/
        
        finalGrid.cell[iRefC[0]].nei.push_back(iRefC[1]);
        finalGrid.cell[iRefC[1]].nei.push_back(iRefC[0]);

        // loop through vertices of face
        for (int v=0; v<f.vtx.size(); ++v)
        {
            // create new point
            ::Point p = pt[ f.vtx[v] ];
            Point pp;
            pp.dim = p.dim;

            int pointIndex;
            //bool ptExists = pointExistsForCreateCells(p, finalGrid.pt, pointIndex);
            bool ptExists = pointExists (pp, fgp, pointIndex);

            if (!ptExists)
            {
                pp.belonging = finalGrid.id;
                finalGrid.pt.push_back (p);
                pointIndex = finalGrid.pt.size() - 1;
                addToPointADT (pp, fgp, pointIndex);
            }

            refFace.vtx.push_back( pointIndex );
        }
                
        for (int i=0; i<neiSize; ++i)
        {
            LRC[i] = NULL;
        }
    }
    
    void findOtherNeiOfGhostsFaces (Grid& finalGrid, PointADT& fgcc)
    {
        for (Face& f: finalGrid.face)
        {
            if (f.bouType == face_t::BOUNDARY)
            {
                if (f.nei.size() == 1)
                {
                    Point cellCent;
                    cellCent.dim = f.cnt;
                    int ic;                    
                    bool pExists = pointExists (cellCent, fgcc, ic);
                    
                    if (pExists)
                    {
                        Cell& cll = finalGrid.cell[ic];
                        
                        f.nei.push_back(ic);
                        finalGrid.cell[f.nei[0]].nei.push_back(ic);
                        cll.nei.push_back(f.nei[0]);
                        int tmp = f.nei[0];
                        f.nei[0] = f.nei[1];
                        f.nei[1] = tmp;

                        cll.face.push_back (&f - &finalGrid.face.front());
                    }
                    else
                    {
                        cout << "could not find in findOtherNeiOfGhostsFaces (...)" << endl;
                        exit(-2);
                    }
                    
                    /*int counter;
                    for (int ic=finalGrid.n_bou_elm; ic<finalGrid.cell.size(); ++ic)
                    {
                        Cell& cll = finalGrid.cell[ic];
                        
                        //if (cll.face.size() != cll.nFaces)
                        //{
                            counter = 0;
                            for (int fv: f.vtx)
                            {
                                for (int cv: cll.vtx)
                                {
                                    if ( cv == fv )
                                    {
                                        ++counter;
                                        break;
                                    }
                                }
                            }

                            if (counter == f.vtx.size())
                            {
                                f.nei.push_back(ic);
                                finalGrid.cell[f.nei[0]].nei.push_back(ic);
                                cll.nei.push_back(f.nei[0]);
                                int tmp = f.nei[0];
                                f.nei[0] = f.nei[1];
                                f.nei[1] = tmp;
                                
                                cll.face.push_back (&f - &finalGrid.face.front());
                                
                                break;
                            }
                        //}
                    }
                    
                    if (counter != f.vtx.size())
                    {
                        cout << "!!! Error: counter != f.vtx.size() in findOtherNeiOfGhostsFaces(...)" << endl;
                        //cout << "face: " << &f - &finalGrid.face.front() << endl;
                        exit(-2);
                    }*/
                }
                else
                {
                    cout << "!!! Error: f.nei.size() != 1 in findOtherNeiOfGhostsFaces(...)" << endl;
                    cout << "!!! Error: f.nei.size() = " << f.nei.size() << endl;
                    exit(-2);
                }
            }
        }
    }
    
    void buildPointADTforFinalGrid (PointADT& pADT, const vector<Grid>& gr)
    {
        pADT.root = new PointADT::Node();
        pADT.root->level = 0;
        pADT.root->p->idx = -1;
        
        for (int d=0; d<ADT_VAR; ++d)
        {
            pADT.root->c[d] = BIG_POS_NUM;
            pADT.root->d[d] = BIG_NEG_NUM;
        }
        
        for (const Grid& g: gr)
        {
            for (const ::Point& p: g.pt)
            {
                for (int d=0; d<ADT_DIM; ++d)
                {
                    pADT.root->c[2*d] = min (p.dim[d], pADT.root->c[2*d]);
                    pADT.root->d[2*d] = max (p.dim[d], pADT.root->d[2*d]);
                }
            }
        }
        
        
        for (int d=0; d<ADT_DIM; ++d)
        {
            pADT.root->c[2*d+1] = pADT.root->c[2*d];
            pADT.root->d[2*d+1] = pADT.root->d[2*d];
        }
    }
    
    void createFinalGrid (Grid& finalGrid, const vector<Grid>& gr, const Grid& newGrid)
    {
        PointADT fgp; // final grid points
        PointADT fgcc; // final grid cell centers
        PointADT ngfc; // new grid face centers
        
        cout << "Building fgp..." << flush;
        buildPointADTforFinalGrid (fgp, gr);
        cout << " done" << endl;
        
        cout << "Building fgcc..." << flush;
        buildPointADTforFinalGrid (fgcc, gr);
        cout << " done" << endl;
        
        
        cout << "Building ngfc..." << flush;
        ngfc.points.reserve (newGrid.face.size());
        
        for (int iF=0; iF<newGrid.face.size(); ++iF)                
        {
            const Face& f = newGrid.face[iF];
            
            ngfc.points.push_back ( ngfc.createADTPoint (f.cnt, f.cnt) );
            ngfc.points.back().idx = iF;
        }
        
        ngfc.build();
        cout << " done" << endl;        
        
        for (const Grid& g: gr)
        {
            for (const Face& f: g.face)
            {
                bool pass = ghostPass (f, g.cell[f.nei[0]].iBlank, g.cell[f.nei[0]].trim);
                if (pass) addGhostCells (f, g.face, g.cell, g.pt, finalGrid, fgp, fgcc);
            }
        }
        
        cout << "added ghost cells of normal grids" << endl;
        //cout << "finalGrid.cell.size() = " << finalGrid.cell.size() << endl;
        //cout << "finalGrid.pt.size() = " << finalGrid.pt.size() << endl;
        
        for (const Face& f: newGrid.face)
        {
            bool pass = ghostPass (f, newGrid.cell[f.nei[0]].iBlank, newGrid.cell[f.nei[0]].trim);
            if (pass) addGhostCells (f, newGrid.face, newGrid.cell, newGrid.pt, finalGrid, fgp, fgcc);
        }
        
        cout << "added ghost cells of new grid" << endl;
        //cout << "finalGrid.cell.size() = " << finalGrid.cell.size() << endl;
        //cout << "finalGrid.pt.size() = " << finalGrid.pt.size() << endl;        
        
        finalGrid.n_bou_elm = finalGrid.cell.size();
        
        for (const Grid& g: gr)
        {
            //const Grid& g = gr[0];
            
            for (const Face& f: g.face)
            {
                bool pass = strictlyInteriorPass (f, g.cell[f.nei[0]].iBlank, g.cell[f.nei[1]].iBlank,
                                                  g.cell[f.nei[0]].trim, g.cell[f.nei[1]].trim);
                if (pass) addCells (f, g.face, g.cell, g.pt, finalGrid, fgp, fgcc);
            }
        }
        
        cout << "added normal grid cells" << endl;
        
        for (const Grid& g: gr)
        {
            //const Grid& g = gr[0];
            
            for (const Face& f: g.face)
            {
                bool pass = intergridPass (f, g.cell[f.nei[0]].iBlank, g.cell[f.nei[1]].iBlank,
                                                  g.cell[f.nei[0]].trim, g.cell[f.nei[1]].trim);
                
                if (pass) addIntergridCells (f, g.face, g.cell, g.pt, finalGrid, newGrid, fgp, fgcc, ngfc);
            }
        }
        
        cout << "added intergrid cells" << endl;
        
        for (const Face& f: newGrid.face)
        {
            if (f.nei.size() == 2)
            {
                bool pass = strictlyInteriorPass (f, newGrid.cell[f.nei[0]].iBlank, newGrid.cell[f.nei[1]].iBlank,
                                                  newGrid.cell[f.nei[0]].trim, newGrid.cell[f.nei[1]].trim);
                
                if (pass) addCells (f, newGrid.face, newGrid.cell, newGrid.pt, finalGrid, fgp, fgcc);
            }
        }
        
        cout << "added new grid cells" << endl;
        
        cout << "Modifying cell Vertices..." << flush;
        modifyCellVertices (finalGrid, newGrid, gr, fgp);
        cout << " done!" << endl;
        
        cout << "Finding Ghosts Faces..." << flush;
        findOtherNeiOfGhostsFaces (finalGrid, fgcc);
        cout << " done!" << endl;
        
        finalGrid.totalNElms = finalGrid.cell.size();
        finalGrid.n_in_elm = finalGrid.totalNElms - finalGrid.n_bou_elm;
        
        cout << "finalGrid.n_bou_elm = " << finalGrid.n_bou_elm << endl;
        cout << "finalGrid.n_in_elm = " << finalGrid.n_in_elm << endl;
    }
}

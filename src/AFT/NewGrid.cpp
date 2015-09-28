#include "AFT.h"

namespace AFT
{
    void createCells (double offsetZ, const vector<Point>& points, Grid& newGrid, const vector<Triangle>& triangles, int phys, int newGridId)
    {
        int iP0, iP1, iP2, iP3;
        vector<Point> points2;
        points2.resize( points.size() );
        CVector na, iv;
        double dp;

        for (unsigned int i=0; i<points.size(); ++i)
        {
            points2[i] = points[i];
            points2[i].dim[2] = offsetZ;
        }
        
        newGrid.n_bou_elm = 2. * triangles.size();        
        newGrid.totalNElms = newGrid.n_bou_elm + triangles.size();
        newGrid.n_in_elm = newGrid.totalNElms - newGrid.n_bou_elm;
        newGrid.cell.resize( newGrid.totalNElms );

        // non-offset boundary elements
        for (unsigned int t=0; t<triangles.size(); ++t)
        {
            newGrid.cell[t].type = elmType_t::TRI;
            newGrid.cell[t].phys = phys;
            newGrid.cell[t].iBlank = iBlank_t::NA;
            newGrid.cell[t].belonging = newGridId;
            newGrid.cell[t].bc = BC::EMPTY;
            newGrid.cell[t].vtx.reserve( static_cast<int>(nVertices_t::TRI) );
            newGrid.cell[t].ghost = true;
            
            iP0 = triangles[t].p[0];
            iP1 = triangles[t].p[1];
            iP2 = triangles[t].p[2];
            
            ::Point p0, p1, p2;
            p0.dim = points[ iP0 ].dim;
            p1.dim = points[ iP1 ].dim;
            p2.dim = points[ iP2 ].dim;
            
            //const Point& p0 = points[ iP0 ];
            //const Point& p1 = points[ iP1 ];
            //const Point& p2 = points[ iP2 ];
            
            int index;
            
            bool p0Exists = pointExistsForCreateCells(p0, newGrid.pt, index);
            
            if (p0Exists)
            {
                newGrid.cell[t].vtx.push_back( index );
            }
            else
            {
                newGrid.pt.push_back( p0 ); newGrid.cell[t].vtx.push_back( newGrid.pt.size()-1 );
            }
            
            bool p1Exists = pointExistsForCreateCells(p1, newGrid.pt, index);
            
            if (p1Exists)
            {
                newGrid.cell[t].vtx.push_back( index );
            }
            else
            {
                newGrid.pt.push_back( p1 ); newGrid.cell[t].vtx.push_back( newGrid.pt.size()-1 );
            }
            
            bool p2Exists = pointExistsForCreateCells(p2, newGrid.pt, index);
            
            if (p2Exists)
            {
                newGrid.cell[t].vtx.push_back( index );
            }
            else
            {
                newGrid.pt.push_back( p2 ); newGrid.cell[t].vtx.push_back( newGrid.pt.size()-1 );
            }

            na = crossP( (p1.dim - p0.dim) , (p2.dim - p0.dim) );
            iv = points2[iP0].dim - p0.dim;
            dp = dotP( na , iv );

            if ( dp < 0. )
            {
                int tmp = newGrid.cell[t].vtx[1];
                newGrid.cell[t].vtx[1] = newGrid.cell[t].vtx[2];
                newGrid.cell[t].vtx[2] = tmp;
                
                /*int tmp = newGrid.face.back().vtx[1];
                newGrid.face.back().vtx[1] = newGrid.face.back().vtx[2];
                newGrid.face.back().vtx[2] = tmp;*/
            }
            else if ( dp == 0. )
            {
                cout << "original equal" << endl;
                cin.ignore();
            }
            
            // assign face
            Face nf;
            nf.vtx.reserve( 3 );
            /*nf.vtx.push_back( newGrid.cell[t].vtx[0] );
            nf.vtx.push_back( newGrid.cell[t].vtx[1] );
            nf.vtx.push_back( newGrid.cell[t].vtx[2] );*/
            nf.vtx.push_back( newGrid.cell[t].vtx[0] );
            nf.vtx.push_back( newGrid.cell[t].vtx[2] );
            nf.vtx.push_back( newGrid.cell[t].vtx[1] );
            nf.bouType = face_t::BOUNDARY;
            newGrid.face.push_back( nf );
            newGrid.cell[t].face.push_back ( newGrid.face.size()-1 );
            newGrid.face.back().nei.push_back (t);
        }

        // offset opposite boundary elements
        for (unsigned int t=triangles.size(); t<newGrid.n_bou_elm; ++t)
        {
            newGrid.cell[t].type = elmType_t::TRI;
            newGrid.cell[t].phys = phys;
            newGrid.cell[t].iBlank = iBlank_t::NA;
            newGrid.cell[t].belonging = newGridId;
            newGrid.cell[t].bc = BC::EMPTY;
            newGrid.cell[t].vtx.reserve( static_cast<int>(nVertices_t::TRI) );
            newGrid.cell[t].ghost = true;
            
            iP0 = triangles[t-triangles.size()].p[0];
            iP1 = triangles[t-triangles.size()].p[1];
            iP2 = triangles[t-triangles.size()].p[2];
            
            ::Point p0, p1, p2;
            p0.dim = points2[ iP0 ].dim;
            p1.dim = points2[ iP1 ].dim;
            p2.dim = points2[ iP2 ].dim;
            
            //const Point& p0 = points2[ iP0 ];
            //const Point& p1 = points2[ iP1 ];
            //const Point& p2 = points2[ iP2 ];
            
            int index;
            
            bool p0Exists = pointExistsForCreateCells(p0, newGrid.pt, index);
            
            if (p0Exists)
            {
                newGrid.cell[t].vtx.push_back( index );
            }
            else
            {
                newGrid.pt.push_back( p0 ); newGrid.cell[t].vtx.push_back( newGrid.pt.size()-1 );
            }
            
            bool p1Exists = pointExistsForCreateCells(p1, newGrid.pt, index);
            
            if (p1Exists)
            {
                newGrid.cell[t].vtx.push_back( index );
            }
            else
            {
                newGrid.pt.push_back( p1 ); newGrid.cell[t].vtx.push_back( newGrid.pt.size()-1 );
            }
            
            bool p2Exists = pointExistsForCreateCells(p2, newGrid.pt, index);
            
            if (p2Exists)
            {
                newGrid.cell[t].vtx.push_back( index );
            }
            else
            {
                newGrid.pt.push_back( p2 ); newGrid.cell[t].vtx.push_back( newGrid.pt.size()-1 );
            }

            na = crossP( (p1.dim - p0.dim) , (p2.dim - p0.dim) );
            iv = points[iP0].dim - p0.dim;
            dp = dotP (na , iv);

            if ( dp > 0. )
            {
                int tmp = newGrid.cell[t].vtx[1];
                newGrid.cell[t].vtx[1] = newGrid.cell[t].vtx[2];
                newGrid.cell[t].vtx[2] = tmp;
                
                /*int tmp = newGrid.face.back().vtx[1];
                newGrid.face.back().vtx[1] = newGrid.face.back().vtx[2];
                newGrid.face.back().vtx[2] = tmp;*/
            }
            else if ( dp == 0. )
            {
                cout << "opposite equal" << endl;
                cin.ignore();
            }

            // assign face
            Face nf;
            nf.vtx.reserve( 3 );
            nf.vtx.push_back( newGrid.cell[t].vtx[0] );
            nf.vtx.push_back( newGrid.cell[t].vtx[1] );
            nf.vtx.push_back( newGrid.cell[t].vtx[2] );
            nf.bouType = face_t::BOUNDARY;
            newGrid.face.push_back( nf );
            newGrid.cell[t].face.push_back ( newGrid.face.size()-1 );
            newGrid.face.back().nei.push_back (t);
        }

        // cells
        for (int t=newGrid.n_bou_elm; t<newGrid.totalNElms; ++t)
        {
            newGrid.cell[t].type = elmType_t::PEN;
            newGrid.cell[t].vtx.reserve( static_cast<int>(nVertices_t::PEN) );
            newGrid.cell[t].face.reserve( static_cast<int>(nFaces_t::PEN) );
            newGrid.cell[t].iBlank = iBlank_t::FIELD;
            newGrid.cell[t].belonging = newGridId;
            newGrid.cell[t].ghost = false;
            
            const Cell& noCell = newGrid.cell[t-newGrid.n_bou_elm];
            const Cell& oCell  = newGrid.cell[t-triangles.size()];

            // assign vertices
            newGrid.cell[t].vtx.push_back( noCell.vtx[0] );
            newGrid.cell[t].vtx.push_back( noCell.vtx[1] );
            newGrid.cell[t].vtx.push_back( noCell.vtx[2] );
            newGrid.cell[t].vtx.push_back( oCell.vtx[0] );
            newGrid.cell[t].vtx.push_back( oCell.vtx[1] );
            newGrid.cell[t].vtx.push_back( oCell.vtx[2] );
            
            // assign face from non-offset boundary element
            newGrid.cell[t].face.push_back ( noCell.face[0] );
            newGrid.face[ noCell.face[0] ].nei.push_back (t);
            
            // assign face from offset boundary element
            newGrid.cell[t].face.push_back ( oCell.face[0] );
            newGrid.face[ oCell.face[0] ].nei.push_back (t);
            
            // internal face 1
            {
                Face nf;
                nf.vtx.reserve( 4 );
                
                nf.vtx.push_back( newGrid.cell[t].vtx[0] );
                nf.vtx.push_back( newGrid.cell[t].vtx[1] );
                nf.vtx.push_back( newGrid.cell[t].vtx[4] );
                nf.vtx.push_back( newGrid.cell[t].vtx[3] );
                
                //nf.vtx.push_back( newGrid.cell[t].vtx[0] );
                //nf.vtx.push_back( newGrid.cell[t].vtx[3] );
                //nf.vtx.push_back( newGrid.cell[t].vtx[4] );
                //nf.vtx.push_back( newGrid.cell[t].vtx[2] );
                //nf.vtx.push_back( newGrid.cell[t].vtx[1] );
                nf.bouType = face_t::INTERIOR;
                int index;
                if (faceExists (nf, newGrid.face, newGrid.pt, index) == false)
                {
                    newGrid.face.push_back( nf );
                    newGrid.cell[t].face.push_back ( newGrid.face.size()-1 );
                    newGrid.face.back().nei.push_back (t);
                }
                else
                {
                    newGrid.cell[t].face.push_back ( index );
                    newGrid.face[ index ].nei.push_back (t);
                }
            }

            /*iP0 = newGrid.face[ newGrid.cell[t].face[2] ].vtx[0];
            iP1 = newGrid.face[ newGrid.cell[t].face[2] ].vtx[1];
            iP2 = newGrid.face[ newGrid.cell[t].face[2] ].vtx[2];
            iP3 = newGrid.face[ newGrid.cell[t].face[2] ].vtx[3];

            na = crossP ( (newGrid.pt[iP1].dim - newGrid.pt[iP0].dim) , (newGrid.pt[iP3].dim - newGrid.pt[iP0].dim) );
            //iv = newGrid.pt[newGrid.cell[t].vtx[1]].dim - newGrid.pt[iP0].dim;
            iv = newGrid.pt[newGrid.cell[t].vtx[2]].dim - newGrid.pt[iP0].dim;
            dp = dotP( na , iv );

            if ( dp > 0. )
            {
                newGrid.face[newGrid.cell[t].face[2]].vtx[1] = iP3;
                newGrid.face[newGrid.cell[t].face[2]].vtx[3] = iP1;
            }
            else if ( dp == 0. )
            {
                cout << "face1 equal" << endl;
                cin.ignore();
            }*/

            // internal face 2
            {
                Face nf;
                nf.vtx.reserve( 4 );
                
                nf.vtx.push_back( newGrid.cell[t].vtx[2] );
                nf.vtx.push_back( newGrid.cell[t].vtx[0] );
                nf.vtx.push_back( newGrid.cell[t].vtx[3] );
                nf.vtx.push_back( newGrid.cell[t].vtx[5] );
                
                //nf.vtx.push_back( newGrid.cell[t].vtx[0] );
                //nf.vtx.push_back( newGrid.cell[t].vtx[3] );
                //nf.vtx.push_back( newGrid.cell[t].vtx[5] );
                //nf.vtx.push_back( newGrid.cell[t].vtx[1] );
                //nf.vtx.push_back( newGrid.cell[t].vtx[2] );
                nf.bouType = face_t::INTERIOR;
                int index;
                if (faceExists (nf, newGrid.face, newGrid.pt, index) == false)
                {
                    newGrid.face.push_back( nf );
                    newGrid.cell[t].face.push_back ( newGrid.face.size()-1 );
                    newGrid.face.back().nei.push_back (t);
                }
                else
                {
                    newGrid.cell[t].face.push_back ( index );
                    newGrid.face[ index ].nei.push_back (t);
                }
            }

            /*iP0 = newGrid.face[newGrid.cell[t].face[3]].vtx[0];
            iP1 = newGrid.face[newGrid.cell[t].face[3]].vtx[1];
            iP2 = newGrid.face[newGrid.cell[t].face[3]].vtx[2];
            iP3 = newGrid.face[newGrid.cell[t].face[3]].vtx[3];

            na = crossP ( (newGrid.pt[iP1].dim - newGrid.pt[iP0].dim) , (newGrid.pt[iP3].dim - newGrid.pt[iP0].dim) );
            //iv = newGrid.pt[newGrid.cell[t].vtx[2]].dim - newGrid.pt[iP0].dim;
            iv = newGrid.pt[newGrid.cell[t].vtx[1]].dim - newGrid.pt[iP0].dim;
            dp = dotP( na , iv );

            if ( dp > 0. )
            {
                newGrid.face[newGrid.cell[t].face[3]].vtx[1] = iP3;
                newGrid.face[newGrid.cell[t].face[3]].vtx[3] = iP1;
            }
            else if ( dp == 0. )
            {
                cout << "face2 equal" << endl;
                cin.ignore();
            }*/

            // internal face 3
            {
                Face nf;
                nf.vtx.reserve( 4 );
                /*nf.vtx.push_back( newGrid.cell[t].vtx[1] );
                nf.vtx.push_back( newGrid.cell[t].vtx[5] );
                nf.vtx.push_back( newGrid.cell[t].vtx[4] );
                nf.vtx.push_back( newGrid.cell[t].vtx[2] );*/
                
                nf.vtx.push_back( newGrid.cell[t].vtx[1] );
                nf.vtx.push_back( newGrid.cell[t].vtx[2] );
                nf.vtx.push_back( newGrid.cell[t].vtx[5] );
                nf.vtx.push_back( newGrid.cell[t].vtx[4] );
                
                nf.bouType = face_t::INTERIOR;
                int index;
                if (faceExists (nf, newGrid.face, newGrid.pt, index) == false)
                {
                    newGrid.face.push_back( nf );
                    newGrid.cell[t].face.push_back ( newGrid.face.size()-1 );
                    newGrid.face.back().nei.push_back (t);
                }
                else
                {
                    newGrid.cell[t].face.push_back ( index );
                    newGrid.face[ index ].nei.push_back (t);
                }
            }

            /*iP0 = newGrid.face[newGrid.cell[t].face[4]].vtx[0];
            iP1 = newGrid.face[newGrid.cell[t].face[4]].vtx[1];
            iP2 = newGrid.face[newGrid.cell[t].face[4]].vtx[2];
            iP3 = newGrid.face[newGrid.cell[t].face[4]].vtx[3];

            na = crossP ( (newGrid.pt[iP1].dim - newGrid.pt[iP0].dim) , (newGrid.pt[iP3].dim - newGrid.pt[iP0].dim) );
            iv = newGrid.pt[newGrid.cell[t].vtx[0]].dim - newGrid.pt[iP0].dim;
            dp = dotP (na , iv );

            if ( dp > 0. )
            {
                newGrid.face[newGrid.cell[t].face[4]].vtx[1] = iP3;
                newGrid.face[newGrid.cell[t].face[4]].vtx[3] = iP1;
            }
            else if ( dp == 0. )
            {
                cout << "face2 equal" << endl;
                cin.ignore();
            }*/
        }
        
        for (Cell& cll: newGrid.cell)
        {
            cll.set_centroid(newGrid.pt);
        }
        
        for (Face& f: newGrid.face)
        {
            f.set_area (newGrid.pt);
            f.set_centroid (newGrid.pt);
            
            if (f.nei.size() == 2)
            {
                newGrid.cell[ f.nei[0] ].nei.push_back( f.nei[1] );
                newGrid.cell[ f.nei[1] ].nei.push_back( f.nei[0] );
            }
            
            if (f.bouType == face_t::BOUNDARY)
            {
                if (f.nei.size() != 2)
                {
                    cout << "f.nei.size() != 2 in createCells(...)" << endl;
                    exit(-2);
                }
                
                int tmp = f.nei[0];
                f.nei[0] = f.nei[1];
                f.nei[1] = tmp;
            }
        }
        
        //for (Face& f: newGrid.face)
        //{
            //if (cll.vol < 0.)
            //{
                //cout << f.area[1] << endl;
              //  exit(-2);
            //}
        //}
        
        //exit(-2);
        
        newGrid.set_elmVolumes();
        
        /*cout << endl;
        for (int ic=newGrid.n_bou_elm; ic<newGrid.cell.size(); ++ic)
        {
            Cell& cll = newGrid.cell[ic];
            
            for (int iF: cll.face)
            {
                Face& f = newGrid.face[iF];
                
                //cout << "cnt = " << f.cnt[1] << endl;
                //cout << "area = " << f.area[1] << endl;
                cout << "tmpd = " << f.cnt[1] * f.area[1] << endl;
            }
            
            cout << "vol = " << cll.vol << endl;
            
            cin.ignore();
        }*/
    }
    
    bool faceExists (const Face& nf, const vector<Face>& face, const vector<::Point>& point, int& index)
    {
        int counter;
        
        for (unsigned int iF=0; iF<face.size(); ++iF)
        {
            const Face& f = face[iF];
            
            counter = 0;
            for (int iV=0; iV<f.vtx.size(); ++iV)
            {
                for (int iNF=0; iNF<nf.vtx.size(); ++iNF)
                {
                    if (f.vtx[iV] == nf.vtx[iNF])
                    {
                        ++counter;
                    }                    
                }
            }
            
            if (counter == f.vtx.size())
            {
                index = iF;
                return true;
            }
        }
        
        return false;
    }
}

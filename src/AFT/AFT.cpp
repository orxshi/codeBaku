#include "AFT.h"

namespace AFT
{
    void aft (vector<Grid>& gr, Grid& finalGrid)
    {
        cout << "Allocating... " << flush;
        int newGridId = 2;
        //int finalGridId = 3;
        double offsetZ = 1.;
        int phys = 1;
        Grid newGrid (gr[0].mainDir, newGridId);        
        //finalGrid.id = finalGridId;
        double aveTriArea;
        //Point meshCenter;
        //meshCenter.dim[0] = 0.5;
        //meshCenter.dim[1] = 0.;
        //meshCenter.dim[2] = 0.;
        vector<Point> points;
        vector<Point> edgeCenters;
        vector<Edge> edges;
        vector<Edge> mesh0Edges;
        vector<Edge> mesh1Edges;
        vector<Triangle> triangles;
        vector<FrontMember> frontList;
        EdgeADT edgeADT;
        EdgeADT edge0ADT;
        EdgeADT edge1ADT;
        EdgeADT edge01ADT;
        TriangleADT triangleADT;
        PointADT pointADT;
        PointADT edgeCenterADT;
        CircleADT circleADT;
        cout << "done!" << endl;
        
        cout << "Trimming/Re-blanking... " << flush;
        gr[0].trimWhoHasFringeNeighbor();
        gr[1].trimWhoHasFringeNeighbor();
        gr[0].trimWhoHasTrimNeighbor (3);
        gr[1].trimWhoHasTrimNeighbor (2);
        cout << "done!" << endl;
        
        cout << "Outputing after trimming/reblanking... " << flush;
        //gr[0].outAllTecplot();
        //gr[1].outAllTecplot();
        gr[0].outAllVTK(0);
        gr[1].outAllVTK(0);
        cout << "done!" << endl;
     
        cout << "Preparing... " << flush;
        setPointsEdges (gr, points, edges, edgeCenters, newGridId);
        
        cout << "edges.size() = " << edges.size() << endl;
        
        createFrontList (edges, frontList, points);
        aveTriArea = getAveTriArea (edges, points);
        cout << "done!" << endl;
        
        cout << "Building trees... " << flush;        
        for (Edge& e: edges)
        {
            if (e.belonging == 0)
            {
                mesh0Edges.push_back (e);
            }
            else if (e.belonging == 1)
            {
                mesh1Edges.push_back (e);
            }
        }
        
        exportToGMSH (points, mesh0Edges, mesh1Edges, gr[0].mainDir); cout << "exported to GMSH" << endl;
        
        edge0ADT.build (points, mesh0Edges); cout << "built edge0ADT" << endl;
        edge1ADT.build (points, mesh1Edges); cout << "built edge1ADT" << endl;
        edgeADT.build (points, edges); cout << "built edgeADT" << endl;
        edge01ADT.build (points, edges); cout << "built edge01ADT" << endl;
        triangleADT.build (edgeADT); cout << "built triangleADT" << endl;
        pointADT.build (points); cout << "built pointADT" << endl;
        edgeCenterADT.build (edgeCenters); cout << "built edgeCenterADT" << endl;
        circleADT.build (edgeADT); cout << "built circleADT" << endl;
        cout << "done!" << endl;
        
        cout << "Advancing front... " << flush;
        advanceFront (frontList, points, aveTriArea, edges, triangles, triangleADT, pointADT, edgeCenterADT, edgeADT, edge01ADT, newGridId, edgeCenters, circleADT);
        cout << "done!" << endl;
        
        cout << "Outputing unflipped triangles... " << flush;        
        outputTrianglesVTK (points, triangles, gr[0].mainDir, "tri.vtk");
        cout << "done!" << endl;
        
        cout << "Neighborhood... " << flush;
        knowParentTriangles (edges, triangles);
        findNeighbors (edges, triangles);
        cout << "done!" << endl;
        
        /*cout << "Flipping triangles... " << flush;
        flip (triangles, edges, points);
        cout << "done!" << endl;
        
        cout << "Flipping triangles... " << flush;
        flip (triangles, edges, points);
        cout << "done!" << endl;*/
        
        cout << "Outputing flipped triangles... " << flush;    
        outputTrianglesVTK (points, triangles, gr[0].mainDir, "triFlip.vtk");
        cout << "done!" << endl;
        
        cout << "Creating cells... " << flush;
        createCells (offsetZ, points, newGrid, triangles, phys, newGridId);
        cout << "done!" << endl;
        
        /*cout << "Outputing new grid... " << flush;
        //newGrid.createOutputDir( gr[0].mainDir );
        //newGrid.outAllTecplot();
        newGrid.outAllVTK(0);
        cout << "done!" << endl;*/
        
        cout << "Creating final grid... " << flush;
        createFinalGrid (finalGrid, gr, newGrid);
        //finalGrid.createOutputDir( gr[0].mainDir );
        cout << "done!" << endl;
        
        /*cout << "Outputing final grid... " << flush;        
        //finalGrid.outAllTecplot();
        finalGrid.outAllVTK(0);
        cout << "done!" << endl;*/
        
        //cout << "Reading input of final grid... " << flush;
        //cout << "FINAL GRID READ INPUT AND PRINT INPUT IS NOT SET IN AFT(...)" << endl;        
        //finalGrid.read_input();
        //finalGrid.printInput();
        //cout << "done!" << endl;
    }
    
    double getAveTriArea (const vector<Edge>& edges, const vector<Point>& points)
    {
        double x,y,size=0.;
        
        for (const Edge& e: edges)
        {
            const Point& t0 = points[ e.t[0] ];
            const Point& t1 = points[ e.t[1] ];
            
            x = t0.dim[0] - t1.dim[0];
            y = t0.dim[1] - t1.dim[1];
            size += pow(x,2) + pow(y,2);
        }
        
        size /= edges.size();
        size *= 2.;
        
        return ( sqrt(3.)/4.*sqrt(size) );
    }
    
    void construct (int iCPX, bool A_CPX_exists, bool B_CPX_exists, int iA_CPX, int iB_CPX, int iA, int iB, vector<FrontMember>& frontList,
             vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, vector<Point>& points, vector<Point>& edgeCenters, PointADT& edgeCenterADT, CircleADT& circleADT)
    {
        int iFrontEdge = frontList.front().edge;
        
        if (!A_CPX_exists)
        {
            Edge tmpEdge = createEdge (iA, iCPX, newGridId, true);
            addToEdgeList (tmpEdge, iA, iCPX, edges, edgeADT, points);
            iA_CPX = edges.size() - 1;
            addToFrontList (iA_CPX, frontList);
            
            Point cntPoint;
            cntPoint.belonging = newGridId;
            cntPoint.dim = 0.5 * (points[tmpEdge.t[0]].dim + points[tmpEdge.t[1]].dim);
            addToPointList (cntPoint, edgeCenters, edgeCenterADT);
        }
        else
        {
            eraseExistingEdgeFromFrontList (iA_CPX, frontList);
        }
        
        if (!B_CPX_exists)
        {
            Edge tmpEdge = createEdge (iB, iCPX, newGridId, true);
            addToEdgeList (tmpEdge, iB, iCPX, edges, edgeADT, points);
            iB_CPX = edges.size() - 1;
            addToFrontList (iB_CPX, frontList);
            
            Point cntPoint;
            cntPoint.belonging = newGridId;
            cntPoint.dim = 0.5 * (points[tmpEdge.t[0]].dim + points[tmpEdge.t[1]].dim);
            addToPointList (cntPoint, edgeCenters, edgeCenterADT);
        }
        else
        {
            eraseExistingEdgeFromFrontList (iB_CPX, frontList);
        }
        
        Triangle tmpTriangle = createTriangle (iFrontEdge, iA_CPX, iB_CPX, edges, points);
        addToTriangleList (triangles, tmpTriangle, triangleADT, points, circleADT, edges);
        
        eraseFromFrontList (frontList);
        sortFrontList (frontList, points, edges);
    }
    
    void exportToGMSH (const vector<Point>& points, const vector<Edge>& mesh0Edges, const vector<Edge>& mesh1Edges, string dir)
    {
        ofstream out;
        dir.append ("/exprtGMSH.geo");
        out.open (dir.c_str());
        
        vector <vector <int> > ptConn ((points.size()+1));
        
        if (out.is_open())
        {
            out << "Mesh.Algorithm = 6;" << endl;
            
            for (int i=0; i<points.size(); ++i)
            {
                out << "Point(";
                out << i;
                out << ") = {";
                out << points[i].dim[0];
                out << ",";
                out << points[i].dim[1];
                out << ",";
                out << points[i].dim[2];
                out << "};" << endl;
            }
            
            //---------------------------------------------
            {
                int start, end;
                
                start = 0;
                end = start + mesh0Edges.size();

                for (int e=start; e<end; ++e) // e for edge
                {
                    int p0 = mesh0Edges[e-start].t[0] ;
                    int p1 = mesh0Edges[e-start].t[1];

                    ptConn[p0].push_back(e);
                    ptConn[p1].push_back(e);
                }

                start = mesh0Edges.size();
                end = start + mesh1Edges.size();

                for (int e=start; e<end; ++e)
                {
                    int p0 = mesh1Edges[e-start].t[0];
                    int p1 = mesh1Edges[e-start].t[1];

                    ptConn[p0].push_back(e);
                    ptConn[p1].push_back(e);
                }
            }
            //---------------------------------------------
            
            out << "Line(";
            out << 0;
            out << ") = {";
            out << mesh0Edges[0].t[0];            
            out << ",";
            out << mesh0Edges[0].t[1];
            out << "};" << endl;
            
            int last = mesh0Edges[0].t[1];
            int iLastEdge = 0;
            
            for (int i=0; i<mesh0Edges.size()-1; ++i)
            {
                // find the edge that last belongs
                for (int j=0; j<2; ++j)
                {
                    if ( ptConn[last][j] != iLastEdge )
                    {
                        // found j for iOtherEdge
                        int iOtherEdge = ptConn[last][j];
                        iLastEdge = iOtherEdge;             
                        
                        // find iOtherEdge.t != it1
                        for (int k=0; k<2; ++k)
                        {
                            if ( mesh0Edges[iOtherEdge].t[k] != last )
                            {                                
                                // found k for iOtherEdge.t != it1                                
                                out << "Line(";
                                out << i+1;
                                out << ") = {";
                                out << last;
                                out << ",";
                                out << mesh0Edges[iOtherEdge].t[k];
                                out << "};" << endl;
                                
                                last = mesh0Edges[iOtherEdge].t[k];                                
                                break;
                            }
                        }
                        
                        break;
                    }
                }
            }
            
            out << "Line(";
            out << mesh0Edges.size();
            out << ") = {";
            out << mesh1Edges[0].t[0];
            out << ",";
            out << mesh1Edges[0].t[1];
            out << "};" << endl;
            
            last = mesh1Edges[0].t[1];
            iLastEdge = mesh0Edges.size();
            
            for (int i=0; i<mesh1Edges.size()-1; ++i)
            {
                // find the edge that last belongs
                for (int j=0; j<2; ++j)
                {
                    if ( ptConn[last][j] != iLastEdge )
                    {
                        // found j for iOtherEdge
                        int iOtherEdge = ptConn[last][j] - mesh0Edges.size();
                        iLastEdge = ptConn[last][j];             
                        
                        // find iOtherEdge.t != it1
                        for (int k=0; k<2; ++k)
                        {
                            if ( mesh1Edges[iOtherEdge].t[k] != last )
                            {
                                // found k for iOtherEdge.t != it1                                
                                out << "Line(";
                                out << i+1 + mesh0Edges.size();
                                out << ") = {";
                                out << last;
                                out << ",";
                                out << mesh1Edges[iOtherEdge].t[k];
                                out << "};" << endl;
                                
                                last = mesh1Edges[iOtherEdge].t[k];
                                break;
                            }
                        }
                        
                        break;
                    }
                }
            }
            
            out << "Line Loop(1) = {";
            for (int i=0; i<mesh0Edges.size(); ++i)
            {
                out << i;
                
                if (i < mesh0Edges.size()-1)
                {
                    out << ",";
                }
            }
            out << "};" << endl;
            
            out << "Line Loop(2) = {";
            for (int i=mesh0Edges.size(); i<mesh0Edges.size()+mesh1Edges.size(); ++i)
            {
                out << i;
                
                if (i < mesh0Edges.size()+mesh1Edges.size()-1)
                {
                    out << ",";
                }
            }
            out << "};" << endl;
            
            out << "Plane Surface(1) = {1,2};";
            
            out.close();
        }
        else
        {
            cout << "could not open file in AFT::exportToGMSH()" << endl;
            exit(-2);
        }
    }
    
    double spacingFnc (double b, double aveTriSize)
    {
        // initially assume uniform mesh
        
        double h = 2. * aveTriSize / b;
        
        return (h/2.);
    }
    
    void Tanemura_Merriam_Helper (int iA, int iB, int& i, vector<Point>& points, deque<int>& pts)
    {
        bool found = false;
        
        CVector center;
        double radius;
        
        Point& A = points[iA];
        Point& B = points[iB];
        
        Point& CPX = points[pts[i]];
        triPtsCircums (CPX.dim, A.dim, B.dim, center, radius);

        cout << "i = " << i << endl;
        cout << "pts[i] = " << pts[i] << endl;
        cout << "CPX.dim.size() = " << CPX.dim.size() << endl;
        cout << "CPX.dim[0] = " << CPX.dim[0] << endl;
        cout << "CPX.dim[1] = " << CPX.dim[1] << endl;
        cout << "radius = " << radius << endl;
                    cout << "center[0] = " << center[0] << endl;
                    cout << "center[1] = " << center[1] << endl;
        
        for (int j=0; j<pts.size(); ++j)
        {
            if (j != i)
            {
                //cout << "points[pts[j]].dim.size() = " << points[pts[j]].dim.size() << endl;
            
                CVector d = points[pts[j]].dim - center;
                
                

                if ( (radius - mag(d)) > 1e-5 )
                {
                    cout << "j = " << j << endl;
                    cout << "pts[j] = " << pts[j] << endl;
                    cout << "points[pts[j]].dim[0] = " << points[pts[j]].dim[0] << endl;
                    cout << "points[pts[j]].dim[1] = " << points[pts[j]].dim[1] << endl;
                    cout << "d[0] = " << d[0] << endl;
                    cout << "d[1] = " << d[1] << endl;
                    cout << "mag(d) = " << mag(d) << endl;
                    cout << "radius2 = " << radius << endl;
                    
                    //cin.ignore();
                    
                    i = j;
                    found = true;
                    break;
                }
            }
        }
                    
        if ( found )
        {
            Tanemura_Merriam_Helper (iA, iB, i, points, pts);
        }
    }
    
    int Tanemura_Merriam (int iA, int iB, vector<Point>& points, deque<int>& pts)
    {
        if (pts.size() == 0)
        {
            cout << "pts.size() == 0 in AFT::Tanemura_Merriam(...)" << endl;
            exit(-2);
        }
        
        cout << "iA = " << iA << endl;
        cout << "iB = " << iB << endl;
        
        int i = 0;
        
        Tanemura_Merriam_Helper (iA, iB, i, points, pts);
        
        return i;
    }
}

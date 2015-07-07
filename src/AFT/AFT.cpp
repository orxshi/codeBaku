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
        double aveCellSize;
        Point meshCenter;
        meshCenter.dim[0] = 0.5;
        meshCenter.dim[1] = 0.;
        meshCenter.dim[2] = 0.;
        vector<Point> points;
        vector<Edge> edges;
        vector<Edge> mesh0Edges;
        vector<Edge> mesh1Edges;
        vector<Triangle> triangles;
        vector<FrontMember> frontList;
        EdgeADT edgeADT;
        EdgeADT edge0ADT;
        EdgeADT edge1ADT;
        TriangleADT triangleADT;
        PointADT pointADT;
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
        setPointsEdges (gr, points, edges, newGridId);
        createFrontList (edges, frontList, points);
        aveCellSize = getAveCellSize (edges, points);
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
        edge0ADT.build (points, mesh0Edges);
        edge1ADT.build (points, mesh1Edges);
        edgeADT.build (points, edges);
        triangleADT.build (edgeADT);
        pointADT.build (points);
        cout << "done!" << endl;
        
        cout << "Advancing front... " << flush;
        advanceFront (frontList, points, aveCellSize, edges, triangles, triangleADT, pointADT, edgeADT, edge0ADT, edge1ADT, newGridId, meshCenter);
        cout << "done!" << endl;
        
        cout << "Neighborhood... " << flush;
        knowParentTriangles (edges, triangles);
        findNeighbors (edges, triangles);
        cout << "done!" << endl;
        
        cout << "Flipping triangles... " << flush;
        flip (triangles, edges, points);
        cout << "done!" << endl;
        
        cout << "Outputing triangles... " << flush;
        outputTriangles (points, triangles);
        cout << "done!" << endl;
        
        cout << "Creating cells... " << flush;
        createCells (offsetZ, points, newGrid, triangles, phys, newGridId);
        cout << "done!" << endl;
        
        cout << "Outputing new grid... " << flush;
        //newGrid.createOutputDir( gr[0].mainDir );
        //newGrid.outAllTecplot();
        newGrid.outAllVTK(0);
        cout << "done!" << endl;
        
        cout << "Creating final grid... " << flush;
        createFinalGrid (finalGrid, gr, newGrid);
        //finalGrid.createOutputDir( gr[0].mainDir );
        cout << "done!" << endl;
        
        cout << "Outputing final grid... " << flush;        
        //finalGrid.outAllTecplot();
        finalGrid.outAllVTK(0);
        cout << "done!" << endl;
        
        cout << "Reading input of final grid... " << flush;
        cout << "FINAL GRID READ INPUT AND PRINT INPUT IS NOT SET IN AFT(...)" << endl;
        exit(-2);
        //finalGrid.read_input();
        //finalGrid.printInput();
        cout << "done!" << endl;
    }
    
    double getAveCellSize (const vector<Edge>& edges, const vector<Point>& points)
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
        
        return size;
    }
}
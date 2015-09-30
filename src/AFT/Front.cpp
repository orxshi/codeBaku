//#include <bits/stl_vector.h>

#include "AFT.h"

namespace AFT
{
    FrontMember::FrontMember()
    {
        newPointChecked = false;
        //cloPtsMaxSize = 30;
    }
    
    void sortFrontList (vector<FrontMember>& frontList, const vector<Point>& points, const vector<Edge>& edges)
    {
        auto cmp = [&] (const FrontMember& fm1, const FrontMember& fm2)
        {
            //cout << "fm2.edge = " << fm2.edge << endl;
            
            const Edge& e1 = edges[ fm1.edge ];
            const Edge& e2 = edges[ fm2.edge ];
            
            const Point& e1t0 = points[ e1.t[0] ];
            const Point& e1t1 = points[ e1.t[1] ];
            const Point& e2t0 = points[ e2.t[0] ];
            const Point& e2t1 = points[ e2.t[1] ];
            
            double mag10 = pow( e1t0.dim[0] - e1t1.dim[0], 2 );
            double mag11 = pow( e1t0.dim[1] - e1t1.dim[1], 2 );
            double mag1  = mag10 + mag11;

            double mag20 = pow( e2t0.dim[0] - e2t1.dim[0], 2 );
            double mag21 = pow( e2t0.dim[1] - e2t1.dim[1], 2 );
            double mag2  = mag20 + mag21;

            if (e1.belonging < e2.belonging)
            {
                return true;
            }
            else if (e1.belonging > e2.belonging)
            {
                return false;
            }
            else if (mag1 < mag2)
            {
                return true;
            }

            return false;
        };
        
        sort(frontList.begin(), frontList.end(), cmp);
        
        /*for (int i=0; i<frontList.size(); ++i)
        {
            cout << edges[frontList[i].edge].belonging << endl;
        }
        
        cin.ignore();*/
    }
    
    void eraseFromFrontList (vector<FrontMember>& frontList)
    {
        if (frontList.size() != 0)
        {
            frontList.erase (frontList.begin());
        }
    }
    
    void addToFrontList (int edge, vector<FrontMember>& frontList)
    {
        FrontMember fm;
        fm.edge = edge;
        
        frontList.push_back(fm);
    }
    
    void eraseExistingEdgeFromFrontList (int ie, vector<FrontMember>& frontList)
    {
        for (unsigned int i=0; i<frontList.size(); ++i)
        {
            if ( frontList[i].edge == ie )
            {
                frontList.erase (frontList.begin() + i);
                return;
            }
        }
    }
    
    void advanceFront (vector<FrontMember>& frontList, vector<Point>& points, double aveTriArea,
                       vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT,
                       PointADT& pointADT, PointADT& edgeCenterADT, EdgeADT& edgeADT, EdgeADT& edge01ADT, int newGridId, vector<Point>& edgeCenters, CircleADT& circleADT)
    {
        #include "advanceFront.h"
        
        double aveEdgeSize = sqrt ( (4./sqrt(3.) * aveTriArea) );
        
        // circumradius of average triangle size
        double rho = triEdgeCircumradius (aveEdgeSize, aveEdgeSize, aveEdgeSize);
        
        while (!frontList.empty())
        {
            FrontMember& frontFirst = frontList.front();
            Edge& frontEdge = edges [ frontFirst.edge ];
            int it0 = frontEdge.t[0];
            int it1 = frontEdge.t[1];
            const Point& t0 = points[ it0 ];
            const Point& t1 = points[ it1 ];
            
            cout << "it0 = " << it0 << endl;
            cout << "it1 = " << it1 << endl;
                        
            // search existing candidate points
            srchCandPts (frontFirst, edges, points, pointADT, candPts, (2.*rho), edgeADT, edge01ADT, triangleADT);
            
            if (candPts.size() == 0)
            {
                cout << "outputing triangles" << endl;
                
                eraseDeadPoints (points, edges, triangles);
                eraseDeadEdges (edges, triangles, points);
                eraseDeadTriangles (triangles, points, edges);
                
                outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
                exit(-2);
            }
            
            // find Delaunay triangle
            int iCandPt = Tanemura_Merriam (it0, it1, points, candPts);
            int iCPX = candPts[iCandPt];
            Point CPX = points[iCPX];
            
            // check whether forming triangle has a circumradius smaller than threshold radius
            bool existingPointPass = checkCircumBound (CPX, t0, t1, rho);
            
            // if circumradius check passed check two forming edges intersections
            if (existingPointPass)
            {
                existingPointPass = checkTwoFormingEdges (CPX, t0, t1, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, edges, edgeADT, points);
            }
            
            // construct triangle if all previous stages are passed
            if (existingPointPass)
            {
                cout << "constructing with existing point" << endl;
                construct (iCPX, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, it0,
                           it1, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points, edgeCenters, edgeCenterADT, circleADT);
            }
            else // if existing point doesn't work
            {
                // get a point normal to front
                Point crP = getNewPt (t0, t1, aveTriArea, edge01ADT, triangleADT);
                
                // no need to check circumradius for uniform grid
                // but will be needed in future
                
                // check two forming edges intersections
                bool newPointPass = checkTwoFormingEdges (crP, t0, t1, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, edges, edgeADT, points);
                
                if (newPointPass)
                {
                    cout << "new point pass in AFT::front(...)" << endl;
                    ptCCInter (crP, circleADT, triangleADT, triangles, frontList, edges, edgeADT, points, pointADT);
                    
                    // add new point to points list
                    addToPointList (crP, points, pointADT);
                    int iCrP = points.size() - 1;
                    
                    cout << "constructing with new point" << endl;
                    // construct triangle if all previous stages are passed                    
                    construct (iCrP, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, it0,
                       it1, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points, edgeCenters, edgeCenterADT, circleADT);
                }
                else
                {
                    cout << "none of the ways work in AFT::advanceFront(...)" << endl;
                
                    cout << "it0 = " << it0 << endl;
                    cout << "it1 = " << it1 << endl;
                    cout << "iCPX = " << iCPX << endl;
                    
                    eraseDeadPoints (points, edges, triangles);
                    eraseDeadEdges (edges, triangles, points);
                    eraseDeadTriangles (triangles, points, edges);

                    outputTrianglesVTK (points, triangles, "../out", "tri.vtk");

                    exit(-2);
                }
            }            
          
            cout << "frontListSize = " << frontList.size() << endl;
        }
    }
}

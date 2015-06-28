#include "AFT.h"

namespace AFT
{
    void F1 (int iCPX, int iCPY, bool& YChecked, bool& XChecked, int terminal, vector<FrontMember>& frontList, vector<Edge>& edges,
             vector<Triangle>& triangles, EdgeADT& edgeADT, TriangleADT& triangleADT, int newGridId, const vector<Point>& points)
    {
        // CPX: Closest point X
        // CPY: Closest point Y
        // t[terminal]   : A
        // t[1-terminal] : B
        bool A_CPX_exists; // t[terminal]-CPX exists

        XChecked = true;
        
        const Point& CPX = points[ iCPX ];
        const Point& CPY = points[ iCPY ];
        FrontMember& frontFirst = frontList.front();
        Edge& frontEdge = edges[ frontFirst.edge ];
        int iA = frontEdge.t[terminal];
        int iB = frontEdge.t[1-terminal];
        const Point& A = points[ iA ];
        const Point& B = points[ iB ];
        
        int iA_CPX = edgeExists (iA, iCPX, points, edges, A_CPX_exists);        

        if (A_CPX_exists)
        {
            Edge& A_CPX = edges[ iA_CPX ];
            F2 (iCPX, iCPY, iA_CPX, YChecked, XChecked, terminal, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points);
        }
        else
        {
            bool intersection = checkEdgeIntersection (A, CPX, edgeADT);

            if (intersection)
            {
                frontFirst.ignore.push_back (iCPX);
                
                if (!YChecked)
                {
                    F1 (iCPY, iCPX, XChecked, YChecked, 1-terminal, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                }
            }
            else
            {
                F3 (iCPX, iCPY, YChecked, XChecked, terminal, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points);
            }
        }
    }
    
    void F2 (int iCPX, int iCPY, int iA_CPX, bool& YChecked, bool& XChecked, int terminal, vector<FrontMember>& frontList,
             vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, const vector<Point>& points)
    {
        // A-CPX already exists so do existingEdge_A_CPX

        bool B_CPX_exists;
        bool tmpBool;
        bool A_B_CPX_exists;
        
        Point CPX = points[ iCPX ];
        Point CPY = points[ iCPY ];        
        int iFrontEdge = frontList.front().edge;
        Edge frontEdge = edges[ iFrontEdge ];
        int iA = frontEdge.t[terminal];
        int iB = frontEdge.t[1-terminal];
        Point A = points[ iA ];
        Point B = points[ iB ];
        Edge A_CPX = edges[ iA_CPX ];

        int iB_CPX = edgeExists (iB, iCPX, points, edges, B_CPX_exists);        

        if (B_CPX_exists)
        {
            Edge B_CPX = edges[ iB_CPX ];
            Triangle tmpTriangle = createTriangle (iFrontEdge, iA_CPX, iB_CPX, edges, points);
            A_B_CPX_exists = triangleIntersect (tmpTriangle, triangleADT, points);

            if (A_B_CPX_exists == false)
            {
                triangles.push_back (tmpTriangle);

                ADT::ADTPoint vec = triangleADT.createADTPoint (tmpTriangle, points);
                vec.idx = triangles.size() - 1;
                triangleADT.insert (vec, triangleADT.root, tmpBool);
                
                eraseExistingEdgeFromFrontList (iA_CPX, frontList);
                eraseExistingEdgeFromFrontList (iB_CPX, frontList);
                eraseFromFrontList (frontList);
            }
            else if (A_B_CPX_exists == true)
            {
                /*
                 * Note that, we don't add a triangle. Consider two cases where front edge is,
                 * 1. An existing edge: Then it would not stay in the front list when
                 *    that triangle formed before.
                 * 2. A newly created edge: Then it should be used for the creation of new triangles
                 *    not for the existing ones.
                */
                
                frontList.front().ignore.push_back ( iCPX );

                if (!YChecked)
                {
                    F1 (iCPY, iCPX, XChecked, YChecked, 1-terminal, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                }
            }
        }
        else // if B_CPX doesn't exist
        {
            bool intersection = checkEdgeIntersection (B, CPX, edgeADT);

            if (intersection)
            {
                frontList.front().ignore.push_back ( iCPX );

                if (!YChecked)
                {
                    F1 (iCPY, iCPX, XChecked, YChecked, 1-terminal, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                }
            }
            else
            {
                if ((CPX.newlyCreated == true) || (B.belonging != CPX.belonging))
                //if ((CPX.belonging == 2) || (frontFirst.t[1-terminal].belonging != CPX.belonging))
                {
                    Edge tmpEdge = createEdge (iB, iCPX, newGridId, true);
                    edges.push_back (tmpEdge);
                    
                    Triangle tmpTriangle = createTriangle (iFrontEdge, iA_CPX, edges.size()-1, edges, points);
                    A_B_CPX_exists = triangleIntersect (tmpTriangle, triangleADT, points);

                    if (!A_B_CPX_exists)
                    {
                        ADT::ADTPoint vec = edgeADT.createADTPoint (points[tmpEdge.t[0]], points[tmpEdge.t[1]]);
                        vec.idx = edges.size() - 1;
                        edgeADT.insert (vec, edgeADT.root, tmpBool);

                        triangles.push_back (tmpTriangle);
                        vec = triangleADT.createADTPoint (tmpTriangle, points);
                        vec.idx = triangles.size() - 1;
                        triangleADT.insert (vec, triangleADT.root, tmpBool);
                        
                        eraseExistingEdgeFromFrontList (iA_CPX, frontList);
                        eraseFromFrontList (frontList);
                        addToFrontList (edges.size()-1, frontList);
                        sortFrontList (frontList, points, edges);
                    }
                    else
                    {
                        edges.pop_back();
                        frontList.front().ignore.push_back ( iCPX );

                        if (!YChecked)
                        {
                            F1 (iCPY, iCPX, XChecked, YChecked, 1-terminal, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                        }
                    }
                }
                else
                {
                    frontList.front().ignore.push_back ( iCPX );

                    if (!YChecked)
                    {
                        F1 (iCPY, iCPX, XChecked, YChecked, 1-terminal, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                    }
                }
            }
        }
    }
    
    void F3 (int iCPX, int iCPY, bool& YChecked, bool& XChecked, int terminal, vector<FrontMember>& frontList, vector<Edge>& edges,
             vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, const vector<Point>& points)
    {
        // A-CPX does not exist

        bool B_CPX_exists;
        bool tmpBool;
        bool A_B_CPX_exists;
        
        const Point CPX = points[ iCPX ];
        const Point CPY = points[ iCPY ];
        int iFrontEdge = frontList.front().edge;
        Edge frontEdge = edges[ iFrontEdge ];
        int iA = frontEdge.t[terminal];
        int iB = frontEdge.t[1-terminal];
        const Point A = points[ iA ];
        const Point B = points[ iB ];
        
        int iB_CPX = edgeExists (iB, iCPX, points, edges, B_CPX_exists);

        if (B_CPX_exists)
        {
            Edge B_CPX = edges[ iB_CPX ];
            Edge tmpEdge1 = createEdge (iA, iCPX, newGridId, true);         
            edges.push_back (tmpEdge1);
            
            Triangle tmpTriangle = createTriangle (iFrontEdge, edges.size()-1, iB_CPX, edges, points);
            A_B_CPX_exists = triangleIntersect (tmpTriangle, triangleADT, points);

            if (!A_B_CPX_exists)
            {                
                ADT::ADTPoint vec = edgeADT.createADTPoint (points[tmpEdge1.t[0]], points[tmpEdge1.t[1]]);
                vec.idx = edges.size() - 1;
                edgeADT.insert (vec, edgeADT.root, tmpBool);

                triangles.push_back(tmpTriangle);
                vec = triangleADT.createADTPoint (tmpTriangle, points);
                vec.idx = triangles.size() - 1;
                triangleADT.insert (vec, triangleADT.root, tmpBool);
                
                eraseExistingEdgeFromFrontList (iB_CPX, frontList);
                eraseFromFrontList (frontList);
                addToFrontList (edges.size()-1, frontList);
                sortFrontList (frontList, points, edges);    
            }
            else
            {
                edges.pop_back();
                
                frontList.front().ignore.push_back ( iCPX );

                if (!YChecked)
                {
                    F1 (iCPY, iCPX, XChecked, YChecked, 1-terminal, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                }
            }
        }
        else
        {
            bool B_CPX_intersects = checkEdgeIntersection (B, CPX, edgeADT);    

            if (B_CPX_intersects)
            {
                frontList.front().ignore.push_back ( iCPX );

                if (!YChecked)
                {
                    F1 (iCPY, iCPX, XChecked, YChecked, 1-terminal, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                }
            }
            else
            {
                Edge A_CPX = createEdge (iA, iCPX, newGridId, true);
                Edge B_CPX = createEdge (iB, iCPX, newGridId, true);
                edges.push_back (A_CPX);
                edges.push_back (B_CPX);

                Triangle tmpTriangle = createTriangle (iFrontEdge, edges.size()-1, edges.size()-2, edges, points);
                A_B_CPX_exists = triangleIntersect (tmpTriangle, triangleADT, points);

                if (!A_B_CPX_exists)
                {
                    ADT::ADTPoint vec = edgeADT.createADTPoint (points[A_CPX.t[0]], points[A_CPX.t[1]]);
                    vec.idx = edges.size() - 2;
                    edgeADT.insert (vec, edgeADT.root, tmpBool);

                    vec = edgeADT.createADTPoint (points[B_CPX.t[0]], points[B_CPX.t[1]]);
                    vec.idx = edges.size() - 1;
                    edgeADT.insert (vec, edgeADT.root, tmpBool);

                    triangles.push_back(tmpTriangle);
                    
                    /*if (triangles.size() == 51)
                                {
                                    cout << "size 51 in f3" << endl;
                                    cout << iA << endl;
                                    cout << iB << endl;
                                    cout << frontEdge.belonging << endl;
                                    cout << iCPX << endl;
                                    cout << CPX.newlyCreated << endl;
                                    cin.ignore();
                                }*/
                    
                    vec = triangleADT.createADTPoint (tmpTriangle, points);
                    vec.idx = triangles.size() - 1;
                    triangleADT.insert (vec, triangleADT.root, tmpBool);

                    eraseFromFrontList (frontList);
                    addToFrontList (edges.size()-2, frontList);
                    addToFrontList (edges.size()-1, frontList);
                    sortFrontList (frontList, points, edges);
                }
                else
                {
                    edges.pop_back();
                    edges.pop_back();

                    frontList.front().ignore.push_back ( iCPX );

                    if (!YChecked)
                    {
                        F1 (iCPY, iCPX, XChecked, YChecked, 1-terminal, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                    }
                }
            }
        }
    }
}

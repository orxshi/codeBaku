/* 
 * File:   AFT.h
 * Author: Orhan Shibliyev
 *
 * Created on August 7, 2014, 11:03 PM
 */

#ifndef AFT_H
#define	AFT_H

#include <algorithm>
#include <memory>
#include <iomanip>
#include "../Grid/Grid.h"

using std::sort;
using std::addressof;
using std::cref;
using std::unique_ptr;
using std::setw;
using std::isfinite;

namespace AFT
{
    struct Edge
    {
        vector<int> t;
        int belonging;
        bool newlyCreated;
        vector<int> nei;
        //int i;

        Edge();
    };
    
    struct Triangle
    {
        vector<int> p;
        vector<int> e;
        vector <int> nei;
        //vector <Edge> neiEdge;
    };
    
    struct FrontMember
    {
        int edge;
        vector<int> ignore;
        int CPfound;
        bool newPointChecked;
        
        FrontMember();
    };
    
    struct EdgeADT : public ADT
    {
        virtual bool compareFunction (const Node* node, const ADTPoint& targetPoint);
        virtual bool doCubesOverlap (const Node* node, const ADTPoint& targetPoint);
        void build (const vector<Point>& points, const vector<Edge>& edges);
        ADTPoint createADTPoint (const Point& a, const Point& b);
    };
    
    struct TriangleADT : public ADT
    {
        virtual bool compareFunction (const Node* node, const ADTPoint& targetPoint);
        virtual bool doCubesOverlap (const Node* node, const ADTPoint& targetPoint);
        void build (const EdgeADT& edgeADT);
        bool localCmpFunc (const ADTPoint& node, const ADTPoint& targetPoint);
        ADTPoint createADTPoint (const Triangle& tri, const vector<Point>& points);
    };
    
    struct PointADT : public ADT
    {
        virtual bool compareFunction (const Node* node, const ADTPoint& targetPoint);
        virtual bool doCubesOverlap (const Node* node, const ADTPoint& targetPoint);        
        void build (const vector<Point>& pointss);
        ADTPoint createADTPoint (const CVector& a, const CVector& b);
    };
    
    // Preset
    void setPointsEdges (const vector<Grid>& gr, vector<Point>& points, vector<Edge>& edges, int newGridId);
    void createFrontList (const vector<Edge>& edges, vector<FrontMember>& frontList, const vector<Point>& points);
    
    // Front    
    void sortFrontList (vector<FrontMember>& frontList, const vector<Point>& points, const vector<Edge>& edges);
    void eraseFromFrontList (vector<FrontMember>& frontList);
    void addToFrontList (int edge, vector<FrontMember>& frontList);
    void eraseExistingEdgeFromFrontList (int ie, vector<FrontMember>& frontList);
    void advanceFront (vector<FrontMember>& frontList, vector<Point>& points, double aveMeshSize,
                       vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT,
                       PointADT& pointADT, EdgeADT& edgeADT, EdgeADT& edge0ADT, EdgeADT& edge1ADT, int newGridId, const Point& meshCenter);
    
    // Edge
    Edge createEdge (int indexA, int indexB, int belonging, bool newlyCreated);
    int edgeExists (const int ip1, const int ip2, const vector<Point>& points, const vector<Edge>& edges, bool& exists);
    bool checkEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT);
    int checkNumberOfEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT);
    void knowParentTriangles (vector<Edge>& edges, const vector<Triangle>& triangles);
    
    // Intersection
    bool doIntersect (const CVector& p1, const CVector& q1, const CVector& p2, const CVector& q2);
    int orientation (const CVector& p, const CVector& q, const CVector& r);
    bool onSegment (const CVector& p, const CVector& q, const CVector& r);
    
    // Triangle
    Triangle createTriangle (int e1, int e2, int e3, const vector<Edge>& edges, const vector<Point>& points);
    bool triangleIntersect (const Triangle& tri, TriangleADT& triangleADT, const vector<Point>& points);
    Triangle createTriangle (int e1, int e2, int e3, const vector<Edge>& edges, const vector<Point>& points);
    double areaTriangle (const Triangle& tri, const vector<Point>& points);
    double charTriangleLength (const Triangle& tri, const vector<Point>& points);
    void findNeighbors (const vector<Edge>& edges, vector<Triangle>& triangles);
    void circleCenter(const CVector& A, const CVector& B, const CVector& C, CVector& cnt, double& radius);
    void flip (vector<Triangle>& triangles, vector<Edge>& edges, const vector<Point>& points);
    void outputTrianglesVTK (const vector<Point>& points, const vector<Triangle>& triangles, string dir);
    void addToTriangleList(vector<Triangle>& triangles, const Triangle& tmpTriangle, TriangleADT& triangleADT, const vector<Point>& points);
    
    // Reblanking
    void fieldToFringe (Grid& gr, Grid& ogr, int crt);
    void fringeToField (Grid& gr, Grid& ogr, int crt);
    
    // Point
    int findClosestPoint (FrontMember& fm, int terminal, const vector<Point>& points, const vector<Edge>& edges, bool& pointFound);
    void createNewPoint (bool& newPointSuccess, vector<FrontMember>& frontList, double aveMeshSize, vector<Point>& points,
                       vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, PointADT& pointADT,
                        EdgeADT& edgeADT, EdgeADT& edge0ADT, EdgeADT& edge1ADT, const Point& meshCenter, const int newGridId);
    void getTwoNormalPoints (int it0, int it1, const vector<Point>& points, Point& crP1, Point& crP2, double pdis);
    bool rayCasting (const Point& farPoint, const Point& p, EdgeADT& edgeADT);
    
    // FModules
    void F1 (int iCPX, int iCPY, bool& YChecked, bool& XChecked, int terminal, vector<FrontMember>& frontList, vector<Edge>& edges,
             vector<Triangle>& triangles, EdgeADT& edgeADT, TriangleADT& triangleADT, int newGridId, const vector<Point>& points);
    void F2 (int iCPX, int iCPY, int iA_CPX, bool& YChecked, bool& XChecked, int terminal, vector<FrontMember>& frontList,
             vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, const vector<Point>& points);
    void F3 (int iCPX, int iCPY, bool& YChecked, bool& XChecked, int terminal, vector<FrontMember>& frontList, vector<Edge>& edges,
             vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, const vector<Point>& points);
    
    double getAveCellSize (const vector<Edge>& edges, const vector<Point>& points);
    void aft (vector<Grid>& gr, Grid& finalGrid);
    void outputTriangles (const vector<Point>& points, const vector<Triangle>& triangles);
    bool faceExists (const Face& nf, const vector<Face>& face, const vector<Point>& point, int& index);
    void createCells (double offsetZ, const vector<Point>& points, Grid& newGrid, const vector<Triangle>& triangles, int phys, int newGridId);
    bool pointExistsForCreateCells(const Point& p, const vector<Point>& points, int& index);
    void createFinalGrid (Grid& finalGrid, const vector<Grid>& gr, const Grid& newGrid);
    void addIntergridCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<Point>& pt, Grid& finalGrid, const Grid& newGrid);
    void addCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<Point>& pt, Grid& finalGrid);
    void modifyCellVertices (Grid& finalGrid, const Grid& newGrid, const vector<Grid>& gr);
};

#endif	/* AFT_H */


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
#include <deque>
#include "../Grid/Grid.h"

using std::sort;
using std::addressof;
using std::cref;
using std::unique_ptr;
using std::setw;
using std::isfinite;
using std::deque;

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
        
        CVector centroid(const vector<Point>& points);
        double qualityScore (const vector<Point>& points, double aveTriArea, bool verbose, bool& passed);
    };
    
    struct FrontMember
    {
        int edge;
        vector<int> ignore;
        int CPfound;
        bool newPointChecked;
        int cloPtsMaxSize;
        deque<int> cloPts;
        
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
        void build ();
        ADTPoint createADTPoint (const CVector& a, const CVector& b);
    };
    
    // Preset
    void setPointsEdges (const vector<Grid>& gr, vector<Point>& points, vector<Edge>& edges, vector<Point>& edgeCenters, int newGridId);
    void createFrontList (const vector<Edge>& edges, vector<FrontMember>& frontList, const vector<Point>& points);
    
    // Front    
    void sortFrontList (vector<FrontMember>& frontList, const vector<Point>& points, const vector<Edge>& edges);
    void eraseFromFrontList (vector<FrontMember>& frontList);
    void addToFrontList (int edge, vector<FrontMember>& frontList);
    void eraseExistingEdgeFromFrontList (int ie, vector<FrontMember>& frontList);
    void advanceFront (vector<FrontMember>& frontList, vector<Point>& points, double aveTriArea,
                       vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT,
                       PointADT& pointADT, PointADT& edgeCenterADT, EdgeADT& edgeADT, EdgeADT& edge01ADT, int newGridId, vector<Point>& edgeCenters);
    
    // Edge
    Edge createEdge (int indexA, int indexB, int belonging, bool newlyCreated);
    int edgeExists (const int ip1, const int ip2, const vector<Point>& points, const vector<Edge>& edges, bool& exists);
    bool checkEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT, const vector<Edge>& edges, const vector<Point>& points, bool& exactMatch, int& result);
    int checkNumberOfEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT);
    void knowParentTriangles (vector<Edge>& edges, const vector<Triangle>& triangles);
    void addToEdgeList (Edge& edge, int iP1, int iP2, vector<Edge>& edges, EdgeADT& edgeADT, const vector<Point>& points);
    
    // Intersection
    bool doIntersect (const CVector& p1, const CVector& q1, const CVector& p2, const CVector& q2, bool& exactMatch);
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
    void outputTrianglesVTK (const vector<Point>& points, const vector<Triangle>& triangles, string dir, string fileName);
    void addToTriangleList(vector<Triangle>& triangles, const Triangle& tmpTriangle, TriangleADT& triangleADT, const vector<Point>& points);
    
    // Reblanking
    void fieldToFringe (Grid& gr, Grid& ogr, int crt);
    void fringeToField (Grid& gr, Grid& ogr, int crt);
    
    // Point
    void findClosestPoint (FrontMember& fm, vector<Edge>& edges, vector<Point>& points);    
    void getTwoNormalPoints (int it0, int it1, const vector<Point>& points, Point& crP1, Point& crP2, double pdis);
    bool rayCasting (const Point& p, EdgeADT& edgeADT);    
    void addToPointList (Point& p, vector<Point>& points, PointADT& pointADT);
    double getPointDistance (double r);
    bool eligible (int iCPX, bool isNewPoint, int iA, int iB, double aveTriArea, double& score, bool& A_CPX_exists, bool& B_CPX_exists, int& iA_CPX, int& iB_CPX, vector<FrontMember>& frontList, vector<Edge>& edges, EdgeADT& edgeADT, EdgeADT& edge01ADT, TriangleADT& triangleADT, vector<Point>& points, PointADT& pointADT, PointADT& edgeCenterADT);
    bool pointsNearby (const CVector& range1, const CVector& range2, PointADT& pointADT, PointADT& edgeCenterADT);
    bool pointExists (const Point& p, PointADT& pointADT, int& result);
    
    // FModules
    //void F1 (int iCPX, int iCPY, bool& YChecked, bool& XChecked, int terminal, vector<FrontMember>& frontList, vector<Edge>& edges,
    //         vector<Triangle>& triangles, EdgeADT& edgeADT, TriangleADT& triangleADT, int newGridId, const vector<Point>& points);
    //void F2 (int iCPX, int iCPY, int iA_CPX, bool& YChecked, bool& XChecked, int terminal, vector<FrontMember>& frontList,
    //         vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, const vector<Point>& points);
    //void F3 (int iCPX, int iCPY, bool& YChecked, bool& XChecked, int terminal, vector<FrontMember>& frontList, vector<Edge>& edges,
    //         vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, const vector<Point>& points);
    
    double getAveTriArea (const vector<Edge>& edges, const vector<Point>& points);
    void aft (vector<Grid>& gr, Grid& finalGrid);
    void outputTriangles (const vector<Point>& points, const vector<Triangle>& triangles);
    bool faceExists (const Face& nf, const vector<Face>& face, const vector<Point>& point, int& index);
    void createCells (double offsetZ, const vector<Point>& points, Grid& newGrid, const vector<Triangle>& triangles, int phys, int newGridId);
    bool pointExistsForCreateCells(const Point& p, const vector<Point>& points, int& index);
    void createFinalGrid (Grid& finalGrid, const vector<Grid>& gr, const Grid& newGrid);
    void addIntergridCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<Point>& pt, Grid& finalGrid, const Grid& newGrid, PointADT& fgp, PointADT& fgcc, PointADT& fgfc);
    void addCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<Point>& pt, Grid& finalGrid, PointADT& fgp, PointADT& fgcc);
    void modifyCellVertices (Grid& finalGrid, const Grid& newGrid, const vector<Grid>& gr);
    void construct (int iCPX, bool isNewPoint, bool A_CPX_exists, bool B_CPX_exists, int iA_CPX, int iB_CPX, int iA, int iB, vector<FrontMember>& frontList,
             vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, const vector<Point>& points, vector<Point>& edgeCenters, PointADT& edgeCenterADT);
};

#endif	/* AFT_H */


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
using std::remove_if;

namespace AFT
{
    struct Point
    {
        bool alive;
        CVector dim;
        int belonging;
        bool newlyCreated;
        vector<int> tri;
        vector<int> e;
        int id;
        
        Point();
        bool eraseParentEdge (int iEdge);
        bool eraseParentTri (int iTri);
    };

    struct Edge
    {
        bool alive;
        vector<int> t;
        int belonging;
        bool newlyCreated;
        vector<int> nei;
        //vector<int> tri;
        int id;

        Edge();
        bool eraseParentTri (int iTri);
    };
    
    struct Triangle
    {   
        bool alive;
        vector<int> p;
        vector<int> e;
        vector <int> nei;
        int id;
        
        Triangle ();
        CVector centroid(const vector<Point>& points);
        double qualityScore (const vector<Point>& points, double aveTriArea, bool verbose, bool& passed);
    };
    
    struct FrontMember
    {
        int edge;
        vector<int> ignore;
        int CPfound;
        bool newPointChecked;
        //int cloPtsMaxSize;
        //deque<int> cloPts;
        
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
        ADTPoint createADTPoint (const Point& p0, const Point& p1, const Point& p2);
    };
    
    struct PointADT : public ADT
    {
        virtual bool compareFunction (const Node* node, const ADTPoint& targetPoint);
        virtual bool doCubesOverlap (const Node* node, const ADTPoint& targetPoint);        
        void build (const vector<Point>& pointss);
        void build ();
        ADTPoint createADTPoint (const CVector& a, const CVector& b);
    };
    
    struct CircleADT : public ADT
    {
        virtual bool compareFunction (const Node* node, const ADTPoint& targetPoint);
        virtual bool doCubesOverlap (const Node* node, const ADTPoint& targetPoint);        
        void build (const EdgeADT& edgeADT);
        ADTPoint createADTPoint (const Triangle& tri, const vector<Point>& points, int idTri);
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
                       PointADT& pointADT, PointADT& edgeCenterADT, EdgeADT& edgeADT, EdgeADT& edge01ADT, int newGridId, vector<Point>& edgeCenters, CircleADT& circleADT);
    
    // Edge
    Edge createEdge (int indexA, int indexB, int belonging, bool newlyCreated);
    int edgeExists (const int ip1, const int ip2, const vector<Point>& points, const vector<Edge>& edges, bool& exists);
    bool checkEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT, const vector<Edge>& edges, const vector<Point>& points, bool& exactMatch, int& result);
    int checkNumberOfEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT);
    void knowParentTriangles (vector<Edge>& edges, const vector<Triangle>& triangles);
    void addToEdgeList (Edge& edge, int iP1, int iP2, vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points);
    void eraseDeadEdges (vector<Edge>& edges, vector<Triangle>& triangles, vector<Point>& points);
    
    // Intersection
    bool doIntersect (const CVector& p1, const CVector& q1, const CVector& p2, const CVector& q2, bool& exactMatch);
    int orientation (const CVector& p, const CVector& q, const CVector& r);
    bool onSegment (const CVector& p, const CVector& q, const CVector& r);
    
    // Triangle
    Triangle createTriangle (int e1, int e2, int e3, const vector<Edge>& edges, const vector<Point>& points);
    bool triangleIntersect (const Triangle& tri, TriangleADT& triangleADT, const vector<Point>& points);
    bool triangleIntersect (const Point& p0, const Point& p1, const Point& p2, TriangleADT& triangleADT);
    Triangle createTriangle (int e1, int e2, int e3, const vector<Edge>& edges, const vector<Point>& points);
    double areaTriangle (const Triangle& tri, const vector<Point>& points);
    double charTriangleLength (const Triangle& tri, const vector<Point>& points);
    void findNeighbors (const vector<Edge>& edges, vector<Triangle>& triangles);
    void triPtsCircums (const CVector& A, const CVector& B, const CVector& C, CVector& cnt, double& radius);
    double triEdgeCircumradius (double a, double b, double c);
    void flip (vector<Triangle>& triangles, vector<Edge>& edges, const vector<Point>& points);
    void outputTrianglesVTK (const vector<Point>& points, const vector<Triangle>& triangles, string dir, string fileName);
    void addToTriangleList(vector<Triangle>& triangles, Triangle& tmpTriangle, TriangleADT& triangleADT, vector<Point>& points, CircleADT& circleADT, vector<Edge>& edges);
    CVector cntTriangle3Pts (const CVector& p0, const CVector& p1, const CVector& p2);
    bool triQuality (const CVector& p0, const CVector& p1, const CVector& p2, double rho);
    void eraseDeadTriangles (vector<Triangle>& triangles, vector<Point>& points, vector<Edge>& edges);
    
    // Reblanking
    void fieldToFringe (Grid& gr, Grid& ogr, int crt);
    void fringeToField (Grid& gr, Grid& ogr, int crt);
    
    // Point
    //void findClosestPoint (FrontMember& fm, vector<Edge>& edges, vector<Point>& points);    
    Point getNewPt (const Point& t0, const Point& t1, double aveTriSize, EdgeADT& edge01ADT, TriangleADT& triangleADT);
    bool rayCasting (const Point& p, EdgeADT& edgeADT);    
    void addToPointList (Point& p, vector<Point>& points, PointADT& pointADT);
    double getPointDistance (double r);
    bool eligible (int iCPX, Point& CPX, bool isNewPoint, int iA, int iB, double aveTriArea, double& score, bool& A_CPX_exists, bool& B_CPX_exists, int& iA_CPX, int& iB_CPX, vector<FrontMember>& frontList, vector<Edge>& edges, EdgeADT& edgeADT, EdgeADT& edge01ADT, TriangleADT& triangleADT, vector<Point>& points, PointADT& pointADT, PointADT& edgeCenterADT, CircleADT& circleADT);
    bool pointsNearby (const CVector& range1, const CVector& range2, PointADT& pointADT, PointADT& edgeCenterADT);
    bool pointExists (const Point& p, PointADT& pointADT, int& result);
    void srchCandPts (FrontMember& fm, vector<Edge>& edges, vector<Point>& points, PointADT& pointADT, deque<int>& candPts, double rho, EdgeADT& edgeADT, EdgeADT& edge01ADT, TriangleADT& triangleADT);
    bool checkTwoFormingEdges (const Point& CPX, const Point& A, const Point& B, bool& A_CPX_exists, bool& B_CPX_exists, int& iA_CPX, int& iB_CPX,
            vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points);
    bool checkCircumBound (const Point& CPX, const Point& A, const Point& B, double rho);    
    void ptCCInter (const Point& CPX, CircleADT& circleADT, TriangleADT& triangleADT, vector<Triangle>& triangles, vector<FrontMember>& frontList, vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points, PointADT& pointADT);
    void eraseDeadPoints (vector<Point>& points, vector<Edge>& edges, vector<Triangle>& triangles);
    
    double getAveTriArea (const vector<Edge>& edges, const vector<Point>& points);
    void aft (vector<Grid>& gr, Grid& finalGrid);
    void outputTriangles (const vector<Point>& points, const vector<Triangle>& triangles);
    bool faceExists (const Face& nf, const vector<Face>& face, const vector<::Point>& point, int& index);
    void createCells (double offsetZ, const vector<Point>& points, Grid& newGrid, const vector<Triangle>& triangles, int phys, int newGridId);
    bool pointExistsForCreateCells(const ::Point& refPoint, const vector<::Point>& points, int& index);
    void createFinalGrid (Grid& finalGrid, const vector<Grid>& gr, const Grid& newGrid);
    void addIntergridCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<::Point>& pt, Grid& finalGrid, const Grid& newGrid, PointADT& fgp, PointADT& fgcc, PointADT& fgfc);
    void addCells (const Face& f, const vector<Face>& face, const vector<Cell>& cell, const vector<::Point>& pt, Grid& finalGrid, PointADT& fgp, PointADT& fgcc);
    void modifyCellVertices (Grid& finalGrid, const Grid& newGrid, const vector<Grid>& gr, PointADT& fgp);
    void construct (int iCPX, bool A_CPX_exists, bool B_CPX_exists, int iA_CPX, int iB_CPX, int iA, int iB, vector<FrontMember>& frontList,
             vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, EdgeADT& edgeADT, int newGridId, vector<Point>& points, vector<Point>& edgeCenters, PointADT& edgeCenterADT, CircleADT& circleADT);
    void exportToGMSH (const vector<Point>& points, const vector<Edge>& mesh0Edges, const vector<Edge>& mesh1Edges, string dir);
    double spacingFnc (double b, double aveTriSize);
    int Tanemura_Merriam (int iA, int iB, vector<Point>& points, deque<int>& pts);
};

#endif	/* AFT_H */


#include "AFT.h"

namespace AFT
{
    void PointADT::build (const vector<Point>& pointss)
    {
        this->points.resize (pointss.size());
        for (unsigned int e=0; e<pointss.size(); ++e)
        {
            this->points[e] = this->createADTPoint (pointss[e].dim, pointss[e].dim);
            this->points[e].idx = e;
        }
        
        ADT::build();
    }
    
    void PointADT::build()
    {
        ADT::build();
    }

    bool PointADT::compareFunction (const Node *node, const ADTPoint& targetPoint)
    {
        return doCubesOverlap (node, targetPoint);
    }

    bool PointADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
    {
        bool insideCube = true;

        for (unsigned int d=0; d<ADT_DIM; ++d)
        {
            if (!(node->p.dim[d*2] <= targetPoint.dim[d*2+1]) || !(node->p.dim[d*2+1] >= targetPoint.dim[d*2]))
            {
                insideCube = false;
                break;
            }
        }

        return insideCube;
    }
    
    ADT::ADTPoint PointADT::createADTPoint (const CVector& a, const CVector& b)
    {
        ADTPoint vec;

        for (unsigned int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = min( a[i], b[i] );
            vec.dim[i*2+1] = max( a[i], b[i] );
        }

        // min = max --> not necessarily

        /** No vertices */

        return vec;
    }
    
    void findClosestPoint (FrontMember& fm, vector<Edge>& edges, vector<Point>& points)
    {
        double dis0;
        double dis1;        
        double tmpDis;
        CVector cVec0;
        CVector cVec1;
        double disI0;
        double disI1;        
        double tmpDisI;
        CVector cVecI0;
        CVector cVecI1;
        
        Edge& frontEdge = edges[fm.edge];
        int it0 = frontEdge.t[0];
        int it1 = frontEdge.t[1];        
        Point& t0= points[it0];
        Point& t1= points[it1];
        
        if (fm.cloPts.size() < fm.cloPtsMaxSize)
        {
            for (unsigned int p=0; p<points.size(); ++p)
            {
                Point& pt = points[p];
                
                if ( p != it1 && p != it0 )
                {
                    cVec0 = pt.dim - t0.dim;
                    cVec1 = pt.dim - t1.dim;
                    dis0 = mag (cVec0);
                    dis1 = mag (cVec1);                    
                    tmpDis = min (dis0, dis1);
                    
                    if ( fm.cloPts.size() == 0 )
                    {
                        fm.cloPts.push_back(p);
                    }
                    else
                    {
                        int size = fm.cloPts.size();
                        
                        for (int i=0; i<fm.cloPts.size(); ++i)
                        {
                            cVecI0 = points[fm.cloPts[i]].dim - t0.dim;
                            cVecI1 = points[fm.cloPts[i]].dim - t1.dim;
                            disI0 = mag (cVecI0);
                            disI1 = mag (cVecI1); 
                            tmpDisI = min (disI0, disI1);

                            if (tmpDis < tmpDisI)
                            {
                                fm.cloPts.insert(fm.cloPts.begin() + i, p);

                                if (fm.cloPts.size() > fm.cloPtsMaxSize)
                                {
                                    fm.cloPts.pop_back();
                                }

                                break;
                            }
                        }
                        
                        if (size == fm.cloPts.size()) // nothing inserted
                        {
                            if (fm.cloPts.size() < fm.cloPtsMaxSize) // still not full
                            {
                                fm.cloPts.push_back(p);
                            }
                        }
                    }
                }
            }
        }
        else if (fm.cloPts.size() == fm.cloPtsMaxSize)
        {
            int p = points.size() - 1;
            Point& pt = points[p];
            
            for (int i: fm.cloPts)
            {
                if (p == i)
                {
                    return;
                }
            }
                
            if (pt.newlyCreated && p != it1 && p != it0 )
            {
                cVec0 = pt.dim - t0.dim;
                cVec1 = pt.dim - t1.dim;
                dis0 = mag (cVec0);
                dis1 = mag (cVec1);                    
                tmpDis = min (dis0, dis1);
                
                for (int i=0; i<fm.cloPts.size(); ++i)
                {
                    cVecI0 = points[fm.cloPts[i]].dim - t0.dim;
                    cVecI1 = points[fm.cloPts[i]].dim - t1.dim;
                    disI0 = mag (cVecI0);
                    disI1 = mag (cVecI1); 
                    tmpDisI = min (disI0, disI1);
                    
                    if (tmpDis < tmpDisI)
                    {
                        fm.cloPts.insert(fm.cloPts.begin() + i, p);
                        fm.cloPts.pop_back();
                        break;
                    }
                }
            }
        }
    }
    
    bool pointsNearby (const CVector& range1, const CVector& range2, PointADT& pointADT, PointADT& edgeCenterADT)
    {
        int result1 = -1;
        int result2 = -1;
        
        ADT::ADTPoint vec = pointADT.createADTPoint (range1, range2);
        
        pointADT.searchForNIntersections = false;
        result1 = pointADT.search (vec);
        
        edgeCenterADT.searchForNIntersections = false;
        result2 = edgeCenterADT.search (vec);

        if (result1 != -1 || result2 != -1)
        {
            return true;
        }
    
        return false;
    }
    
    bool pointExists (const Point& p, PointADT& pointADT, int& result)
    {
        result = -1;
        
        CVector meshDis;
        meshDis[0] = 1e-10;
        meshDis[1] = 1e-10;
        meshDis[2] = 1e-10;
        
        ADT::ADTPoint vec = pointADT.createADTPoint (p.dim-meshDis, p.dim+meshDis);
        
        /*cout << "p.dim[0] = " << p.dim[0] << endl;
        cout << "p.dim[1] = " << p.dim[1] << endl;
        cout << "p.dim[2] = " << p.dim[2] << endl;
        
        cout << "vec.dim[0] = " << vec.dim[0] << endl;
        cout << "vec.dim[1] = " << vec.dim[1] << endl;
        cout << "vec.dim[2] = " << vec.dim[2] << endl;*/
        
        pointADT.searchForNIntersections = false;
        result = pointADT.search (vec);

        if (result != -1)
        {
            return true;
        }
    
        return false;
    }
    
    double getPointDistance (double r)
    {
        return (0.5 * r * 1.73205080757); // weird number is tan(60)
    }
    
    

    void getTwoNormalPoints (int it0, int it1, const vector<Point>& points, Point& crP1, Point& crP2, double pdis)
    {
        double dx = points[it1].dim[0] - points[it0].dim[0];
        double dy = points[it1].dim[1] - points[it0].dim[1];
        double dz = 0.;
        double length = sqrt( pow(dx,2) + pow(dy,2) + pow(dz,2) );

        CVector normal1;
        normal1[0] = -dy;
        normal1[1] = dx;
        normal1[2] = 0.;

        CVector normal2;
        normal2[0] = dy;
        normal2[1] = -dx;
        normal2[2] = 0.;

        normal1 = norm (normal1);
        normal2 = norm (normal2);

        CVector center;
        center[0] = 0.5 * (points[it0].dim[0] + points[it1].dim[0]);
        center[1] = 0.5 * (points[it0].dim[1] + points[it1].dim[1]);
        center[2] = 0.;

        //Point crP1;
        //Point crP2;

        pdis = getPointDistance (length);

        crP1.dim = center + (pdis * normal1);
        crP2.dim = center + (pdis * normal2);
    }
    
    bool rayCasting (const Point& p, EdgeADT& edgeADT)
    {
        // determines whether point, p is inside domain of edgeADT or not
        // returns true if inside
        // returns false if outside
        
        Point farPoint;
        farPoint.dim[0] = edgeADT.root->d[0] + fabs(edgeADT.root->d[0]); // 2*xmax
        farPoint.dim[1] = edgeADT.root->d[1] + fabs(edgeADT.root->d[1]); // 2*ymax
        farPoint.dim[2] = 0.;
        
        int nInter = checkNumberOfEdgeIntersection (farPoint, p, edgeADT);

        if (nInter % 2 == 0) // even (outside)
        {
            return false;
        }
        else // odd (inside)
        {
            return true;
        }
    }
    
    void addToPointList (Point& p, vector<Point>& points, PointADT& pointADT)
    {
        bool tempBool;
        
        points.push_back(p);
        ADT::ADTPoint vec = pointADT.createADTPoint (p.dim, p.dim);
        vec.idx = points.size() - 1;
        pointADT.insert (vec, pointADT.root, tempBool);
    }
    
    bool eligible (int iCPX, bool isNewPoint, int iA, int iB, double aveTriArea, double& score, bool& A_CPX_exists, bool& B_CPX_exists, int& iA_CPX, int& iB_CPX, vector<FrontMember>& frontList, vector<Edge>& edges, EdgeADT& edgeADT, EdgeADT& edge01ADT, TriangleADT& triangleADT, vector<Point>& points, PointADT& pointADT, PointADT& edgeCenterADT)
    {
        // don't forget to pop back edges
        
        bool passed;
        score = BIG_POS_NUM;
        
        bool A_CPX_intersects;
        bool B_CPX_intersects;
        
        const Point& CPX = points[ iCPX ];
        FrontMember& frontFirst = frontList.front();
        int iFrontEdge = frontFirst.edge;
        //Edge& frontEdge = edges[ iFrontEdge ];
        //int iA = frontEdge.t[terminal];
        //int iB = frontEdge.t[1-terminal];
        const Point& A = points[ iA ];
        const Point& B = points[ iB ];
        
        A_CPX_intersects = checkEdgeIntersection (A, CPX, edgeADT, edges, points, A_CPX_exists, iA_CPX);
        //cout << "A_CPX_intersects = " << A_CPX_intersects << endl;
        if (A_CPX_intersects)
        {
            /*cout << "A_CPX_exists = " << A_CPX_exists << endl;
            cout << "edges[iA_CPX].t[0] = " << edges[iA_CPX].t[0] << endl;
            cout << "edges[iA_CPX].t[0] = " << edges[iA_CPX].t[1] << endl;
            cout << "iA_CPX = " << iA_CPX << endl;*/
        }
        if (A_CPX_intersects && A_CPX_exists || !A_CPX_intersects)
        {
            B_CPX_intersects = checkEdgeIntersection (B, CPX, edgeADT, edges, points, B_CPX_exists, iB_CPX);
            //cout << "B_CPX_intersects = " << B_CPX_intersects << endl;
            if (B_CPX_intersects)
            {
                /*cout << "B_CPX_exists = " << B_CPX_exists << endl;
                cout << "edges[iB_CPX].t[0] = " << edges[iB_CPX].t[0] << endl;
                cout << "edges[iB_CPX].t[0] = " << edges[iB_CPX].t[1] << endl;
                cout << "iB_CPX = " << iB_CPX << endl;*/
            }
            if (B_CPX_intersects && B_CPX_exists || !B_CPX_intersects)
            {
                Edge A_CPX;
                Edge B_CPX;
                
                if (A_CPX_exists) { A_CPX = edges[ iA_CPX ]; }
                if (B_CPX_exists) { B_CPX = edges[ iB_CPX ]; }
                
                if (!A_CPX_exists)
                {
                    int dummyID = -1;
                    A_CPX = createEdge (iA, iCPX, dummyID, true);
                    edges.push_back (A_CPX);
                    iA_CPX = edges.size() - 1;
                    
                    
                }
                if (!B_CPX_exists)
                {
                    int dummyID = -1;
                    B_CPX = createEdge (iB, iCPX, dummyID, true);
                    edges.push_back (B_CPX);
                    iB_CPX = edges.size() - 1;
                }
                
                Triangle tmpTriangle = createTriangle (iFrontEdge, iA_CPX, iB_CPX, edges, points);
                bool A_B_CPX_intersects = triangleIntersect (tmpTriangle, triangleADT, points);
                Point tmpCntPnt;
                tmpCntPnt.dim = tmpTriangle.centroid (points);
                bool centInsideDomain = rayCasting (tmpCntPnt, edge01ADT);  
                
                if (!A_B_CPX_intersects && centInsideDomain)
                {
                    if (isNewPoint)
                    {
                        if ( rayCasting (CPX, edge01ADT) )
                        {
                            CVector meshDis;
                            meshDis[0] = 1e-10;
                            meshDis[1] = 1e-10;
                            meshDis[2] = 0.;
                            
                            int iExistPoint;
                            if ( !(pointExists (CPX, pointADT, iExistPoint)) )
                            //if ( !(pointExists (CPX.dim-meshDis, CPX.dim+meshDis, pointADT)) )
                            {
                                double edgeAveTri = sqrt ( (4./sqrt(3.) * aveTriArea) );
                                double pdisAveTri = getPointDistance (edgeAveTri);
                                double pdis = getPointDistance ( mag (A.dim - B.dim) );
                                meshDis[0] = 0.5 * pdis;
                                meshDis[1] = 0.5 * pdis;
                                meshDis[2] = 0.;
                                
                                /*if (iA == 60 && iB == 57)
                                {
                                    cout << "here here" << endl;
                                    cout << pointsNearby (CPX.dim-meshDis, CPX.dim+meshDis, pointADT, edgeCenterADT) << endl;
                                    exit(-2);
                                }*/
                                                                
                                if ( !(pointsNearby (CPX.dim-meshDis, CPX.dim+meshDis, pointADT, edgeCenterADT)) )
                                {
                                    score = tmpTriangle.qualityScore (points, aveTriArea, false, passed);
                                    if (!A_CPX_exists) {edges.pop_back();}
                                    if (!B_CPX_exists) {edges.pop_back();}
                                    
                                    return passed;
                                }
                            }
                        }
                    }
                    else
                    {
                        score = tmpTriangle.qualityScore (points, aveTriArea, false, passed);
                        if (!A_CPX_exists) {edges.pop_back();}
                        if (!B_CPX_exists) {edges.pop_back();}
                        
                        return passed;
                    }
                }
                
                if (!A_CPX_exists) {edges.pop_back();}
                if (!B_CPX_exists) {edges.pop_back();}
            }            
        }                
        
        return false;
    }
}


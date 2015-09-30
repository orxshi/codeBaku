#include "AFT.h"

namespace AFT
{
    Point::Point ()
    {
        belonging = -1;
        newlyCreated = false;
        alive = true;
    }

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
            if (!(node->p->dim[d*2] <= targetPoint.dim[d*2+1]) || !(node->p->dim[d*2+1] >= targetPoint.dim[d*2]))
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
    
    void moveNewPt (Point& p, const vector<int>& iPts, double srchRegion, const vector<Point>& points)
    {
        CVector d, Dif;
        double radius;
        double dif;
        
        radius = srchRegion;
        
        for (int i=0; i<iPts.size(); ++i)
        {
            d = p.dim - points[iPts[i]].dim;
            dif = radius - mag(d);
            Dif = dif * norm(d);
            p.dim += Dif;
        }
    }
    
    bool srchNearbyPts (const Point& p, const vector<Point>& points, int t0, int t1, PointADT& pointADT, double srchRegion, vector<int>& iPts)
    {
        //int it0 = edges[fm.edge].t[0];
        //int it1 = edges[fm.edge].t[1];        
        //Point& t0 = points[it0];
        //Point& t1 = points[it1];        
        //double d = mag (t0.dim - t1.dim);
        //double srchRegion = d * coef;
        
        CVector range1;
        range1[0] = p.dim[0] - srchRegion;
        range1[1] = p.dim[1] - srchRegion;
        range1[2] = 0.;
        
        CVector range2;
        range2[0] = p.dim[0] + srchRegion;
        range2[1] = p.dim[1] + srchRegion;
        range2[2] = 0.;
        
        ADT::ADTPoint vec = pointADT.createADTPoint (range1, range2);
        
        pointADT.searchForNIntersections = true;
        pointADT.search (vec);
        iPts.clear();
        //iPts = pointADT.ids;
        
        CVector d;
        double radius;
        double dif;
        
        radius = srchRegion;
        
        cout << "pointADT.ids.size() = " << pointADT.ids.size()  << endl;
        
        for (int i=0; i<pointADT.ids.size(); ++i)
        {
            if (pointADT.ids[i] != t0 && pointADT.ids[i] != t1)
            {
                cout << "pointADT.ids[i] != t0 && pointADT.ids[i] != t1" << endl;
                cout << "pointADT.ids[i] = " << pointADT.ids[i] << endl;
                
                d = p.dim - points[pointADT.ids[i]].dim;
                dif = radius - mag(d);
                
                cout << "p.dim[1] = " << p.dim[0] << endl;
                cout << "p.dim[2] = " << p.dim[1] << endl;
                
                cout << "points[pointADT.ids[i]].dim[0] = " << points[pointADT.ids[i]].dim[0] << endl;
                cout << "points[pointADT.ids[i]].dim[1] = " << points[pointADT.ids[i]].dim[1] << endl;
                
                cout << "radius = " << radius << endl;
                cout << "mag(d) = " << mag(d) << endl;
                cout << "dif = " << dif << endl;

                if (dif > 0.)
                {
                    // inside circle
                    cout << "inside circle" << endl;
                    iPts.push_back (pointADT.ids[i]);
                }
            }
        }
        
        
        cout << "iPts.size() = " << iPts.size()  << endl;
        
        if (iPts.size() > 0)
        {
            return true;
        }
        else
        {
            return false;
        }
        
        
        
        /*if (result == -1)
        {
            return false;
        }
        else
        {
            return true;
        }*/
        
        /*candPts.clear();
        candPts.push_back (pointADT.ids[0]);
        
        for (int i=1; i<pointADT.ids.size(); ++i)
        {
            Point p = points [pointADT.ids[i]];
            
            double d1 = mag (p.dim - t0.dim);
            double d2 = mag (p.dim - t1.dim);
            if ( min(d1,d2) < candPts.back() )
            {
                candPts.push_back (pointADT.ids[i]);
            }
            else
            {
                candPts.push_front (pointADT.ids[i]);
            }
        }*/
    }
    
    void srchCandPts (FrontMember& fm, vector<Edge>& edges, vector<Point>& points, PointADT& pointADT, deque<int>& candPts, double rho, EdgeADT& edgeADT, EdgeADT& edge01ADT, TriangleADT& triangleADT)
    {
        int it0 = edges[fm.edge].t[0];
        int it1 = edges[fm.edge].t[1];        
        Point& t0 = points[it0];
        Point& t1 = points[it1];        
        double d = mag (t0.dim - t1.dim);
        double srchRegion = rho;
        
        CVector range1;
        range1[0] = min (t0.dim[0], t1.dim[0]);
        range1[1] = min (t0.dim[1], t1.dim[1]);
        range1[0] -= srchRegion;
        range1[1] -= srchRegion;
        range1[2] = 0.;
        
        CVector range2;
        range2[0] = max (t0.dim[0], t1.dim[0]);
        range2[1] = max (t0.dim[1], t1.dim[1]);
        range2[0] += srchRegion;
        range2[1] += srchRegion;
        range2[2] = 0.;
        
        ADT::ADTPoint vec = pointADT.createADTPoint (range1, range2);
        
        pointADT.searchForNIntersections = true;
        pointADT.search (vec);
        
        if (pointADT.ids.size() == 0)
        {
            cout << "pointADT.ids.size() == 0 in AFT::srchCandPts(...)" << endl;
            //exit(-2);
        }
        
        candPts.clear();
        
        for (int i=0; i<pointADT.ids.size(); ++i)
        {
            const Point& p = points [pointADT.ids[i]];
            
            if (pointADT.ids[i] != it1 && pointADT.ids[i] != it0 )
            {
                Point tmpCntPnt;
                tmpCntPnt.dim = cntTriangle3Pts (t0.dim, t1.dim, p.dim);
                bool cntInsideDomain = rayCasting (tmpCntPnt, edge01ADT);
                
                if (cntInsideDomain)
                {
                    bool t0_P_exist;
                    int dummyInt;
                    
                    bool t0_p_inter = checkEdgeIntersection (t0, p, edgeADT, edges, points, t0_P_exist, dummyInt);
                    
                    if (!t0_p_inter || t0_P_exist)
                    {
                        bool t1_P_exist;
                        bool t1_p_inter = checkEdgeIntersection (t1, p, edgeADT, edges, points, t1_P_exist, dummyInt);
                        
                        if (!t1_p_inter || t1_P_exist)
                        {
                            if (t0_P_exist && t1_P_exist)
                            {
                                // check if triangle exists
                                if (triangleIntersect (p, t0, t1, triangleADT))
                                {
                                    continue;
                                }
                            }
                        
                            double d1 = mag (p.dim - t0.dim);
                            double d2 = mag (p.dim - t1.dim);

                            if ( candPts.empty() )
                            {
                                candPts.push_back (pointADT.ids[i]);
                            }
                            else
                            {
                                if ( min(d1,d2) < candPts.back() )
                                {
                                    candPts.push_back (pointADT.ids[i]);
                                }
                                else
                                {
                                    candPts.push_front (pointADT.ids[i]);
                                }
                            }
                        }
                        else
                        {
                            cout << "BBBpointADT.ids.size() = 0 in AFT::srchCandPts(...)" << endl;
                            
                        }
                    }
                    else
                    {
                        cout << "AAApointADT.ids.size() = 0 in AFT::srchCandPts(...)" << endl;
                        
                    }
                }
                else
                {
                    cout << "CCCpointADT.ids.size() = 0 in AFT::srchCandPts(...)" << endl;
                }
            }
        }
        
        if (candPts.size() == 0)
        {
            cout << "no cand points found in AFT::srchCandPts(...)" << endl;
            //exit(-2);
        }
    }
    
    /*void findClosestPoint (FrontMember& fm, vector<Edge>& edges, vector<Point>& points)
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
    }*/
    
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
    
    

    Point getNewPt (const Point& t0, const Point& t1, double aveTriSize, EdgeADT& edge01ADT, TriangleADT& triangleADT)
    {
        Point tmpCntPnt;
        bool cntInsideDomain;
        Point crP;
        
        double dx = t1.dim[0] - t0.dim[0];
        double dy = t1.dim[1] - t0.dim[1];
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
        center[0] = 0.5 * (t0.dim[0] + t1.dim[0]);
        center[1] = 0.5 * (t0.dim[1] + t1.dim[1]);
        center[2] = 0.;
        
        double s = spacingFnc (length, aveTriSize);

        crP.dim = center + (s * normal1);
        
        tmpCntPnt.dim = cntTriangle3Pts (t0.dim, t1.dim, crP.dim);
        cntInsideDomain = rayCasting (tmpCntPnt, edge01ADT);
        
        if (cntInsideDomain)
        {
            if (triangleIntersect (crP, t0, t1, triangleADT))
            {
                crP.dim = center + (s * normal2);
                
                tmpCntPnt.dim = cntTriangle3Pts (t0.dim, t1.dim, crP.dim);
                cntInsideDomain = rayCasting (tmpCntPnt, edge01ADT);
                
                if (cntInsideDomain)
                {
                    return crP;
                }
                else
                {
                    cout << "cannot get a proper new point in AFT::getNewPt(...)" << endl;
                    exit(-2);
                }
            }
            else
            {
                return crP;
            }
        }
        else
        {
            crP.dim = center + (s * normal2);
            
            tmpCntPnt.dim = cntTriangle3Pts (t0.dim, t1.dim, crP.dim);
            cntInsideDomain = rayCasting (tmpCntPnt, edge01ADT);
            
            if (cntInsideDomain)
            {
                return crP;
            }
            else
            {
                cout << "cannot get a proper new point in AFT::getNewPt(...)" << endl;
                exit(-2);
            }
        }
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
    
    /*bool eligible (int iCPX, Point& CPX, bool isNewPoint, int iA, int iB, double aveTriArea, double& score, bool& A_CPX_exists, bool& B_CPX_exists, int& iA_CPX, int& iB_CPX, vector<FrontMember>& frontList, vector<Edge>& edges, EdgeADT& edgeADT, EdgeADT& edge01ADT, TriangleADT& triangleADT, vector<Point>& points, PointADT& pointADT, PointADT& edgeCenterADT, CircleADT& circleADT)
    {
        // don't forget to pop back edges
        
        bool passed;
        score = BIG_POS_NUM;
        
        bool A_CPX_intersects;
        bool B_CPX_intersects;
        
        //Point& CPX = points[ iCPX ];
        FrontMember& frontFirst = frontList.front();
        int iFrontEdge = frontFirst.edge;
        //Edge& frontEdge = edges[ iFrontEdge ];
        //int iA = frontEdge.t[terminal];
        //int iB = frontEdge.t[1-terminal];
        const Point& A = points[ iA ];
        const Point& B = points[ iB ];
        
        A_CPX_intersects = checkEdgeIntersection (A, CPX, edgeADT, edges, points, A_CPX_exists, iA_CPX);
        
        if (A_CPX_intersects && A_CPX_exists || !A_CPX_intersects)
        {
            B_CPX_intersects = checkEdgeIntersection (B, CPX, edgeADT, edges, points, B_CPX_exists, iB_CPX);
            
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
                                //double pdis = getPointDistance ( mag (A.dim - B.dim) );
                                double pdis = mag (A.dim - B.dim);
                                
                                cout << "iA = " << iA << endl;
                                cout << "iB = " << iB << endl;
                                
                                meshDis[0] = 1.0 * pdis;
                                meshDis[1] = 1.0 * pdis;
                                //meshDis[0] = pdisAveTri;
                                //meshDis[1] = pdisAveTri;
                                meshDis[2] = 0.;
                                
                                //
                                ADT::ADTPoint vecC;
                                
                                vecC.dim[0] = CPX.dim[0];
                                vecC.dim[2] = CPX.dim[1];
                                vecC.dim[4] = CPX.dim[2];
                                
                                vecC.dim[1] = vecC.dim[0];
                                vecC.dim[3] = vecC.dim[2];
                                vecC.dim[5] = vecC.dim[4];
                                
                                circleADT.searchForNIntersections = false;
                                int res = circleADT.search (vecC);
                                if (res == -1)
                                {
                                    score = tmpTriangle.qualityScore (points, aveTriArea, false, passed);
                                    if (!A_CPX_exists) {edges.pop_back();}
                                    if (!B_CPX_exists) {edges.pop_back();}
                                    cout << "new point passed" << endl;
                                    cout << "score = " << score << endl;
                                    cout << "CPX.dim[0]; = " << CPX.dim[0] << endl;
                                    cout << "CPX.dim[1]; = " << CPX.dim[1] << endl;
                                    return passed;
                                }
                                //
                                
                                /*vector<int> iPts;
                                if ( !srchNearbyPts (CPX, points, iA, iB, pointADT, pdis, iPts) )
                                //if ( !(pointsNearby (CPX.dim-meshDis, CPX.dim+meshDis, pointADT, edgeCenterADT)) )
                                {
                                    cout << "no points nearby" << endl;
                                    
                                    score = tmpTriangle.qualityScore (points, aveTriArea, false, passed);
                                    cout << "score = " << score << endl;
                                    cout << "passed = " << passed << endl;
                                    if (!A_CPX_exists) {edges.pop_back();}
                                    if (!B_CPX_exists) {edges.pop_back();}
                                    
                                    return passed;
                                }
                                else
                                {
                                    cout << "points nearby" << endl;
                                    moveNewPt (CPX, iPts, pdis, points);
                                    
                                    if ( !srchNearbyPts (CPX, points, iA, iB, pointADT, pdis, iPts) )
                                    {
                                        if (rayCasting (CPX, edge01ADT))
                                        {
                                            int dummyID2 = -1;
                                            
                                            Edge A_CPX2 = createEdge (iA, iCPX, dummyID2, true);
                                            edges.push_back (A_CPX2);
                                            int iA_CPX2 = edges.size() - 1;
                                            
                                            Edge B_CPX2 = createEdge (iB, iCPX, dummyID2, true);
                                            edges.push_back (B_CPX2);
                                            int iB_CPX2 = edges.size() - 1;
                                            
                                            Triangle tmpTriangle2 = createTriangle (iFrontEdge, iA_CPX2, iB_CPX2, edges, points);
                                            
                                            score = tmpTriangle2.qualityScore (points, aveTriArea, false, passed);
                                            
                                            edges.pop_back();
                                            edges.pop_back();
                                            if (!A_CPX_exists) {edges.pop_back();}
                                            if (!B_CPX_exists) {edges.pop_back();}
                                            
                                            return passed;
                                        }
                                    }
                                }*/
                            /*}
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
    }*/
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*bool eligible (int iCPX, Point& CPX, bool isNewPoint, int iA, int iB, double aveTriArea, double& score, bool& A_CPX_exists, bool& B_CPX_exists, int& iA_CPX, int& iB_CPX, vector<FrontMember>& frontList, vector<Edge>& edges, EdgeADT& edgeADT, EdgeADT& edge01ADT, TriangleADT& triangleADT, vector<Point>& points, PointADT& pointADT, PointADT& edgeCenterADT, CircleADT& circleADT)
    {
        bool passed;
        
        bool A_CPX_intersects;
        bool B_CPX_intersects;
        
        FrontMember& frontFirst = frontList.front();
        int iFrontEdge = frontFirst.edge;
        const Point& A = points[ iA ];
        const Point& B = points[ iB ];
        
        //--check circumbound-----------------------------------------------------
        double aveEdgeSize = sqrt ( (4./sqrt(3.) * aveTriArea) );
        bool inCircumBound = triQuality (A.dim, B.dim, CPX.dim, aveEdgeSize/2.)
        
        if (!inCircumBound)
        {
            return false;
        }
        
        //--check A_CPX-----------------------------------------------------
        A_CPX_intersects = checkEdgeIntersection (A, CPX, edgeADT, edges, points, A_CPX_exists, iA_CPX);
        
        if (A_CPX_intersects)
        {
            if (!A_CPX_exists)
            {
                return false;
            }
        }
        
        //--check B_CPX-----------------------------------------------------
        B_CPX_intersects = checkEdgeIntersection (B, CPX, edgeADT, edges, points, B_CPX_exists, iB_CPX);
        
        if (B_CPX_intersects)
        {
            if (!B_CPX_exists)
            {
                return false;
            }
        }
            
        if (isNewPoint)
        {
            if ( rayCasting (CPX, edge01ADT) )
            {
                    // check if new point intersects a circumcircle
                    ADT::ADTPoint vecC;

                    vecC.dim[0] = CPX.dim[0];
                    vecC.dim[2] = CPX.dim[1];
                    vecC.dim[4] = CPX.dim[2];

                    vecC.dim[1] = vecC.dim[0];
                    vecC.dim[3] = vecC.dim[2];
                    vecC.dim[5] = vecC.dim[4];

                    circleADT.searchForNIntersections = false;
                    int res = circleADT.search (vecC);
                    if (res == -1)
                    {
                        score = tmpTriangle.qualityScore (points, aveTriArea, false, passed);
                        if (!A_CPX_exists) {edges.pop_back();}
                        if (!B_CPX_exists) {edges.pop_back();}
                        cout << "new point passed" << endl;
                        cout << "score = " << score << endl;
                        cout << "CPX.dim[0]; = " << CPX.dim[0] << endl;
                        cout << "CPX.dim[1]; = " << CPX.dim[1] << endl;
                        return passed;
                    }
            }
        }
                       
                      
        
        return true;
    }*/
    
    void ptCCInter (const Point& CPX, CircleADT& circleADT, TriangleADT& triangleADT, vector<Triangle>& triangles, vector<FrontMember>& frontList, vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points, PointADT& pointADT)    
    {
        // checks if new point intersects circumcircles
        // delete corresponding triangles from triangles list and tree
        
        ADT::ADTPoint vecC;

        vecC.dim[0] = CPX.dim[0];
        vecC.dim[2] = CPX.dim[1];
        vecC.dim[4] = CPX.dim[2];

        vecC.dim[1] = vecC.dim[0];
        vecC.dim[3] = vecC.dim[2];
        vecC.dim[5] = vecC.dim[4];

        cout << "searching in circleADT in ptCCInter(...)" << endl;
        
        circleADT.searchForNIntersections = true;
        circleADT.search (vecC);
        
        cout << "circleADT.ids.size() = " << circleADT.ids.size() << endl;
        cout << "done searching in circleADT in ptCCInter(...)" << endl;
        
        if (circleADT.ids.size() != 0)
        {
            cout << "circleADT.ids.size() != 0" << endl;
            
            
        
            // new point intersected at least one circumcircle            
            // circleADT will return ids of triangles
            
            // pick a triangle
            for (int t: circleADT.ids)
            {
                if (t == -1)
                {
                    cout << "circleADT.ids[] = -1 which should not be in AFT::ptccInter(...)" << endl;
                    exit (-2);
                }
            
                cout << "circleADT.ids = " << t << endl;
            
                // check neighbors via edges
                for (int e: triangles[t].e)
                {
                    cout << "e = " << e << endl;
                    
                    if (edges[e].alive)
                    {                
                        // if triangle has no neighbor on this edge
                        if (edges[e].nei.size() == 1)
                        {
                            cout << "triangle has no neighbor on this edge" << endl;
                        
                            // kill the edge
                            edges[e].alive = false;
                            // remove edge from edgeADT
                            cout << "removing from edgeADT" << endl;
                            bool success = edgeADT.removeViaID (e);
                            cout << "done removing from edgeADT" << endl;
                            if (!success)
                            {
                                cout << "could not remove from edgeADT 1 in AFT::ptCCInter(...)" << endl;
                                cout << "e = " << e << endl;
                                cout << "t = " << t << endl;
                                cout << "triangles[t].e[0] = " << triangles[t].e[0] << endl;
                                cout << "triangles[t].e[1] = " << triangles[t].e[1] << endl;
                                cout << "triangles[t].e[2] = " << triangles[t].e[2] << endl;
                                cout << "triangles[t].p[0] = " << triangles[t].p[0] << endl;
                                cout << "triangles[t].p[1] = " << triangles[t].p[1] << endl;
                                cout << "triangles[t].p[2] = " << triangles[t].p[2] << endl;
                                cout << "edges[e].t[0] = " << edges[e].t[0] << endl;
                                cout << "edges[e].t[1] = " << edges[e].t[1] << endl;
                                exit (-2);
                            }
                            // remove the edge from front list
                            for (int fr=0; fr<frontList.size(); ++fr)
                            {
                                if (frontList[fr].edge == e)
                                {
                                    frontList.erase (frontList.begin() + fr);
                                    break;
                                }
                            }
                            
                            cout << "notify correspoding points about the killing of edge" << endl;
                            // notify correspoding points about the killing of edge
                            
                            cout << "edges[e].t[0] = " << edges[e].t[0] << endl;
                            cout << "edges[e].t[1] = " << edges[e].t[1] << endl;
                            cout << "points.size() = " << points.size() << endl;
                            
                            points[edges[e].t[0]].eraseParentEdge (e);
                            points[edges[e].t[1]].eraseParentEdge (e);
                            
                            cout << "done notify correspoding points about the killing of edge" << endl;
                        }
                        else if (edges[e].nei.size() == 2)
                        {
                            cout << "triangle has a neighbor on this edge" << endl;
                        
                            // find neighbor
                            for (int n: edges[e].nei)
                            {                           
                                if (n != t)
                                {
                                    bool removeNeiToo = false;
                                
                                    for (int tt: circleADT.ids)
                                    {
                                        if (tt == -1)
                                        {
                                            cout << "circleADT.ids[] = -1 which should not be in AFT::ptccInter(...)" << endl;
                                            exit (-2);
                                        }
                                    
                                        // if triangle has a neighbor which is also going to be removed
                                        if (tt == n)
                                        {
                                            removeNeiToo = true;
                                            // kill the edge
                                            edges[e].alive = false;
                                            // remove edge from edgeADT
                                            bool success = edgeADT.removeViaID (e);
                                            
                                            /*if (e == 217)
                                            {
                                                cout << "deleteing edge 217" << endl;
                                                cout << "edges[e].t[0] = " << edges[e].t[0] << endl;
                                                cout << "edges[e].t[1] = " << edges[e].t[1] << endl;
                                                exit(-2);
                                            }*/
                                            
                                            if (!success) {cout << "could not remove from edgeADT 2 in AFT::ptCCInter(...)" << endl; exit (-2);}
                                            // notify correspoding points about the killing of edge
                                            points[edges[e].t[0]].eraseParentEdge (e);
                                            points[edges[e].t[1]].eraseParentEdge (e);
                                            
                                            break;
                                        }
                                    }
                                    
                                    // update neighbour of edge                                        
                                    for (int nn=0; nn<edges[e].nei.size(); ++nn)
                                    {
                                        if (edges[e].nei[nn] == t)
                                        {
                                            edges[e].nei.erase (edges[e].nei.begin() + nn);
                                            break;
                                        }
                                    }
                                    
                                    if (edges[e].nei.size() >= 2)
                                    {
                                        cout << "edges[e].nei.size() >= 2 in AFT::ptCCInter(...)" << endl;
                                        cout << "edges[e].nei.size() = " << edges[e].nei.size() << endl;
                                        exit (-2);
                                    }
                                    
                                    // if triangle has a neighbor which is NOT going to be removed
                                    if (!removeNeiToo)
                                    {
                                        // add edge to front list
                                        addToFrontList (e, frontList);
                                    }
                                    
                                    break;
                                }
                            }
                        }
                        else
                        {
                            cout << "undefined situation in AFT::ptCCInter(...)" << endl;
                            exit (-2);
                        }
                    }
                }
                
                cout << "killing triangle" << endl;
                
                // kill the triangle
                triangles[t].alive = false;
                // remove triangle from triangleADT
                bool success = triangleADT.removeViaID (t);
                if (!success) {cout << "could not remove from triangleADT in AFT::ptCCInter(...)" << endl; exit (-2);}
                // notify correspoding points about the killing of triangle
                points[triangles[t].p[0]].eraseParentTri (t);
                points[triangles[t].p[1]].eraseParentTri (t);
                points[triangles[t].p[2]].eraseParentTri (t);
                // notify correspoding edges about the killing of triangle
                edges[triangles[t].e[0]].eraseParentTri (t);
                edges[triangles[t].e[1]].eraseParentTri (t);
                edges[triangles[t].e[2]].eraseParentTri (t);
                // remove circle from circleADT
                success = circleADT.removeViaID (t);
                if (!success) {cout << "could not remove from circleADT in AFT::ptCCInter(...)" << endl; exit (-2);}
            }
        }
        
        // kill unwanted points and remove them from trees
        for (int p=0; p<points.size(); ++p)
        {
            if (points[p].e.empty() && points[p].alive == true)
            {
                points[p].alive = false;
                cout << "removing from pointADT" << endl;                
                
                bool success = pointADT.removeViaID (p);
                
                cout << "done removing from pointADT" << endl;
                if (!success) {cout << "could not remove from pointADT in AFT::ptCCInter(...)" << endl; exit (-2);}
            }
        }
        
        cout << "point intersects circumcircle(s) in AFT::ptCCInter(...)" << endl;
        //exit(-2);
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    bool checkTwoFormingEdges (const Point& CPX, const Point& A, const Point& B, bool& A_CPX_exists, bool& B_CPX_exists, int& iA_CPX, int& iB_CPX,
            vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points)
    {
        //--check A_CPX-----------------------------------------------------
        bool A_CPX_intersects = checkEdgeIntersection (A, CPX, edgeADT, edges, points, A_CPX_exists, iA_CPX);
        
        if (A_CPX_intersects)
        {
            if (!A_CPX_exists)
            {
                cout << "A_CPX_intersects in AFT::checkTwoFormingEdges(...)" << endl;
                return false;
            }
        }
        
        //--check B_CPX-----------------------------------------------------
        bool B_CPX_intersects = checkEdgeIntersection (B, CPX, edgeADT, edges, points, B_CPX_exists, iB_CPX);
        
        if (B_CPX_intersects)
        {
            if (!B_CPX_exists)
            {
                cout << "B_CPX_intersects in AFT::checkTwoFormingEdges(...)" << endl;
                cout << "CPX.dim[0] = " << CPX.dim[0] << endl;
                cout << "CPX.dim[1] = " << CPX.dim[1] << endl;
                return false;
            }
        }
        
        return true;
    }
    
    
    bool checkCircumBound (const Point& CPX, const Point& A, const Point& B, double rho)
    {
        bool inCircumBound = triQuality (A.dim, B.dim, CPX.dim, rho);
        
        if (!inCircumBound)
        {
            cout << "!inCircumBound in AFT::checkCircumBound(...)" << endl;
            cout << "rho = " << rho << endl;
            return false;
        }
        
        return true;
    }
    
    bool Point::eraseParentEdge (int iEdge)
    {
        for (int ie=0; ie<e.size(); ++ie)
        {
            int ee = e[ie];
        
            if (ee == iEdge)
            {
                e.erase (e.begin() + ie);
                return true;
            }
        }
        
        return false;
    }
    
    bool Point::eraseParentTri (int iTri)
    {
        for (int it=0; it<tri.size(); ++it)        
        {
            int t = tri[it];
        
            if (t == iTri)
            {
                tri.erase (tri.begin() + it);
                return true;
            }
        }
        
        return false;
    }
    
    void eraseDeadPoints (vector<Point>& points, vector<Edge>& edges, vector<Triangle>& triangles)
    {
        for (int ip=0; ip<points.size(); ++ip)
        {
            Point& p = points[ip];
            
            p.id = ip;
        }
        
        points.erase (remove_if (points.begin(), points.end(), [](Point& p) { return p.alive == false; }), points.end());
        
        /*for (int ip=0; ip<points.size(); ++ip)
        {
            Point& p = points[ip];
        
            if (p.alive == false)
            {
                points.erase(points.begin() + ip);
            }
        }*/
        
        for (int ip=0; ip<points.size(); ++ip)
        {
            Point& p = points[ip];
            
            for (int e: points[p.id].e)
            {            
                Edge& edge = edges[e];
                
                for (int ee=0; ee<edge.t.size(); ++ee)
                {
                    if (edge.t[ee] == p.id)
                    {
                        edge.t[ee] = ip;
                        break;
                    }
                }
            }
            
            for (int t: points[p.id].tri)
            {            
                Triangle& tri = triangles[t];
                
                for (int pp=0; pp<tri.p.size(); ++pp)
                {
                    if (tri.p[pp] == p.id)
                    {
                        tri.p[pp] = ip;
                        break;
                    }
                }
            }
        }
        
        for (int ip=0; ip<points.size(); ++ip)
        {
            Point& p = points[ip];
            
            p.id = ip;
        }
    }
}


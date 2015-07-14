#include "AFT.h"

namespace AFT
{
    void PointADT::build (const vector<Point>& pointss)
    {
        this->points.resize (pointss.size());
        for (unsigned int e=0; e<pointss.size(); ++e)
        {
            this->points[e] = this->createADTPoint (pointss[e].dim, pointss[e].dim);
        }
        
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
    
    int findClosestPoint (FrontMember& fm, int terminal, const vector<Point>& points, const vector<Edge>& edges, bool& pointFound)
    {
        pointFound = false;
        double dis = BIG_POS_NUM;
        double tmpDis;
        bool cancel;
        int index = -1;
        bool cond;
        
        const Edge& frontEdge = edges[ fm.edge ];
        
        const int it0 = frontEdge.t[terminal];
        const int it1 = frontEdge.t[1-terminal];
        
        const Point& t0= points[ it0 ];
        const Point& t1= points[ it1 ];
        
        auto func = [&] (unsigned int p)
        {
            cancel = false;
            for (const int& igr: fm.ignore)
            {
                if (p == igr)
                {
                    cancel = true;
                    break;
                }
            }

            if (cancel == true)
            {
                //continue;
                return;
            }
            else
            {
                tmpDis = sqrt( pow(t0.dim[0] - points[p].dim[0],2) + pow(t0.dim[1] - points[p].dim[1],2) );

                if ( tmpDis < dis )
                {
                    dis = tmpDis;
                    index = p;
                }
            }
        };
        
        /*if (frontEdge.newlyCreated)
        {
            for (unsigned int p=0; p<points.size(); ++p)
            {
                if ( p != it1 && p != it0 )
                {
                    func(p);
                }
            }
        }
        else
        {
            for (int p=0; p<points.size(); ++p)
            {
                if (points[p].belonging != frontEdge.belonging)
                {
                    func(p);
                }
            }
        }*/
        
        for (unsigned int p=0; p<points.size(); ++p)
        {
            if ( p != it1 && p != it0 )
            {
                func(p);
            }
        }
        
        if (index != -1)
        {
            ++fm.CPfound;
            pointFound = true;
        }
        
        /*if ( t0.newlyCreated )
        {
            for (int p=0; p<points.size(); ++p)
            {
                if (points[p].belonging != t0.belonging)
                {
                    if ( p != it1 )
                    {
                        cancel = false;
                        for (const int& igr: fm.ignore)
                        {
                            if ( p == igr )
                            {
                                cancel = true;
                                break;
                            }
                        }
                        
                        if (cancel == true)
                        {
                            continue;
                        }
                        else
                        {
                            tmpDis = sqrt( pow(t0.dim[0] - points[p].dim[0],2) + pow(t0.dim[1] - points[p].dim[1],2) );

                            if ( tmpDis < dis )
                            {
                                dis = tmpDis;
                                index = p;
                            }
                            
                        }
                    }
                }
            }
            
            if (index != -1)
            {
                ++fm.CPfound;
                pointFound = true;
            }
        }
        else
        {
            for (int p=0; p<points.size(); ++p)
            {
                if ( (p != it1) && (p != it0) )
                {
                    cancel = false;
                    for (const int& igr: fm.ignore)
                    {
                        if ( p == igr )
                        {
                            cancel = true;
                            break;
                        }
                    }
                    
                    if (cancel == true)
                    {
                        continue;
                    }
                    else
                    {
                        tmpDis = sqrt( pow(t0.dim[0] - points[p].dim[0],2) + pow(t0.dim[1] - points[p].dim[1],2) );

                        if ( tmpDis < dis )
                        {
                            dis = tmpDis;
                            index = p;
                        }
                    }
                }
            }
            
            if (index != -1)
            {
                ++fm.CPfound;
                pointFound = true;
            }
        }*/
        
        return index;
    }
    
    bool pointsNearby (const CVector& range1, const CVector& range2, PointADT& pointADT)
    {
        int result = -1;
        
        ADT::ADTPoint vec = pointADT.createADTPoint (range1, range2);

        pointADT.searchForNIntersections = false;
        result = pointADT.search (vec);

        if (result != -1)
        {
            return true;
        }
    
        return false;
    }
    
    bool pointExists (const CVector& p, const vector<Point>& points)
    {
        for (unsigned int i=0; i<points.size(); ++i)
        {
            if ( fabs(p[0] - points[i].dim[0]) < 0.0001 )
            {
                if ( fabs(p[1] - points[i].dim[1]) < 0.0001 )
                {
                    return true;
                }
            }
        }

        return false;
    }
    
    double getPointDistance (double r)
    {
        return (0.5 * r * 1.73205080757); // weird number is tan(60)
    }
    
    void createTriWithNewPoint (Point& crP, vector<FrontMember>& frontList, vector<Point>& points,
                        vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, PointADT& pointADT,
                        EdgeADT& edgeADT, const int newGridId)
    {
        FrontMember& frontFirst = frontList.front();
        int iFrontEdge = frontFirst.edge;
        Edge& frontEdge = edges[iFrontEdge];
        int it0 = frontEdge.t[0];
        int it1 = frontEdge.t[1];
        
        crP.belonging = newGridId;
        crP.newlyCreated = true;            
        addToPointList (crP, points, pointADT);
        int iCrP = points.size() - 1;

        Edge tmpEdge1 = createEdge (iCrP, it0, newGridId, true);
        addToEdgeList (tmpEdge1, iCrP, it0, edges, edgeADT, points);
        int iTmpEdge1 = edges.size() - 1;

        Edge tmpEdge2 = createEdge (iCrP, it1, newGridId, true);
        addToEdgeList (tmpEdge2, iCrP, it1, edges, edgeADT, points);
        int iTmpEdge2 = edges.size() - 1;

        Triangle tmpTriangle;
        tmpTriangle = createTriangle (iFrontEdge, iTmpEdge1, iTmpEdge2, edges, points);
        addToTriangleList (triangles, tmpTriangle, triangleADT, points);

        eraseFromFrontList (frontList);
        addToFrontList (iTmpEdge2, frontList);
        addToFrontList (iTmpEdge1, frontList);
        sortFrontList (frontList, points, edges);

        return;
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
    
    bool candidateNewPoint (Point& crP, double& scoreNP, double aveTriArea, vector<FrontMember>& frontList, double pdisAveTri, vector<Point>& points,
                        vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, PointADT& pointADT,
                        EdgeADT& edgeADT, EdgeADT& edge01ADT)
    {
        Triangle tmpTriangle;
        Point crP1;
        Point crP2;
        double pdis;
        
        FrontMember& frontFirst = frontList.front();
        int iFrontEdge = frontFirst.edge;
        Edge& frontEdge = edges[iFrontEdge];
        int it0 = frontEdge.t[0];
        int it1 = frontEdge.t[1];
        
        bool newPointSuccess = true;

        getTwoNormalPoints (it0, it1, points, crP1, crP2, pdis);

        points.push_back (crP1); int iCrP1 = points.size() - 1;
        points.push_back (crP2); int iCrP2 = points.size() - 1;

        if (frontEdge.newlyCreated)
        {
            tmpTriangle.p.push_back (it0);
            tmpTriangle.p.push_back (it1);
            tmpTriangle.p.push_back (iCrP1);

            bool triIntersects = triangleIntersect (tmpTriangle, triangleADT, points);

            if (!triIntersects)
            {
                bool crP1InsideDomain = rayCasting (crP1, edge01ADT);
                if (crP1InsideDomain)
                {
                    crP = crP1;
                }
                else
                {
                    newPointSuccess = false;
                }
            }
            else
            {
                bool crP2InsideDomain = rayCasting (crP2, edge01ADT);
                if (crP2InsideDomain)
                {
                    crP = crP2;
                }
                else
                {
                    newPointSuccess = false;
                }
            }
        }
        else
        {
            bool crP1InsideDomain = rayCasting (crP1, edge01ADT);
            bool crP2InsideDomain = rayCasting (crP2, edge01ADT);

            if (crP1InsideDomain)
            {
                crP = crP1;
            }
            else if (crP2InsideDomain)
            {
                crP = crP2;
            }
            else
            {
                newPointSuccess = false;
            }
        }

        points.pop_back(); // erase crP1
        points.pop_back(); // erase crP2

        if (newPointSuccess)
        {
            newPointSuccess = false;

            bool doesPointExist = pointExists (crP.dim, points);

            if (!doesPointExist)
            {
                CVector meshDis;
                meshDis[0] = 0.7 * pdisAveTri;
                meshDis[1] = 0.7 * pdisAveTri;
                meshDis[2] = 0.;

                bool pointNearby = pointsNearby (crP.dim-meshDis, crP.dim+meshDis, pointADT);

                if (!pointNearby)
                {
                    bool intersection = checkEdgeIntersection (crP, points[it0], edgeADT);

                    if (!intersection)
                    {
                        intersection = checkEdgeIntersection (crP, points[it1], edgeADT);

                        if (!intersection)
                        {
                            newPointSuccess = true;
                            points.push_back(crP);
                            int iNP = points.size() - 1;

                            Triangle tmpTriangle;

                            tmpTriangle.p.push_back (it0);
                            tmpTriangle.p.push_back (it1);
                            tmpTriangle.p.push_back (iNP);

                            scoreNP = tmpTriangle.qualityScore (points, aveTriArea);
                            
                            points.pop_back();
                        }
                    }
                }
            }
        }
        
        return newPointSuccess;
    }
    
    void addToPointList (Point& p, vector<Point>& points, PointADT& pointADT)
    {
        bool tempBool;
        
        points.push_back(p);
        ADT::ADTPoint vec = pointADT.createADTPoint (p.dim, p.dim);
        vec.idx = points.size() - 1;
        pointADT.insert (vec, pointADT.root, tempBool);
    }
    
    
}


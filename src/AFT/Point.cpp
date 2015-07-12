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
        
        if (frontEdge.newlyCreated)
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
    
    void createNewPoint (bool& newPointSuccess, vector<FrontMember>& frontList, double aveMeshSize, vector<Point>& points,
                        vector<Edge>& edges, vector<Triangle>& triangles, TriangleADT& triangleADT, PointADT& pointADT,
                        EdgeADT& edgeADT, EdgeADT& edge0ADT, EdgeADT& edge1ADT, const Point& meshCenter, const int newGridId)
    {
        Triangle tmpTriangle;
        bool tempBool;
        ADT::ADTPoint vec;
        
        FrontMember& frontFirst = frontList.front();
        int iFrontEdge = frontFirst.edge;
        Edge& frontEdge = edges[iFrontEdge];
        int it0 = frontEdge.t[0];
        int it1 = frontEdge.t[1];
        Point farPoint;
        farPoint.dim[0] = pointADT.root->d[0] + fabs(pointADT.root->d[0]); // 2*xmax
        farPoint.dim[1] = pointADT.root->d[1] + fabs(pointADT.root->d[1]); // 2*ymax
        farPoint.dim[2] = 0.;        
        
        if ( frontFirst.newPointChecked == false )
        {
            /*double dx = points[it1].dim[0] - points[it0].dim[0];
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

            Point crP;
            Point crP1;
            Point crP2;
            
            double pdis = getPointDistance (length);
            
            crP1.dim = center + (pdis * normal1);
            crP2.dim = center + (pdis * normal2);  */
            
            Point crP;
            Point crP1;
            Point crP2;
            double pdis;
            
            getTwoNormalPoints (it0, it1, points, crP1, crP2, pdis);
                        
            points.push_back (crP1);
            points.push_back (crP2);
            
            int iCrP1 = points.size() - 2;
            int iCrP2 = points.size() - 1;

            if (frontEdge.newlyCreated)
            //if (frontEdge.belonging == 2)
            {
                tmpTriangle.p.push_back (it0);
                tmpTriangle.p.push_back (it1);
                tmpTriangle.p.push_back (iCrP1);
                
                bool triIntersects = triangleIntersect (tmpTriangle, triangleADT, points);
                
                if (!triIntersects)
                {
                    //crP = crP1;
                    bool crP1InsideDomain = rayCasting (farPoint, crP1, edgeADT);
                    if (crP1InsideDomain)
                    {
                        crP = crP1;
                    }
                    else
                    {
                        cout << "crP1 is outside domain in useNewPoints(...)" << endl;
                        frontFirst.newPointChecked = true;
                        newPointSuccess = false;
                        points.pop_back(); // erase crP1
                        points.pop_back(); // erase crP2
                        return;
                    }
                }
                else
                {
                    //crP = crP2;
                    bool crP2InsideDomain = rayCasting (farPoint, crP2, edgeADT);
                    if (crP2InsideDomain)
                    {
                        crP = crP2;
                    }
                    else
                    {
                        cout << "crP2 is outside domain in useNewPoints(...)" << endl;
                        frontFirst.newPointChecked = true;
                        newPointSuccess = false;
                        points.pop_back(); // erase crP1
                        points.pop_back(); // erase crP2
                        return;
                    }
                }
            }
            else
            {
                /*bool p1e1 = checkEdgeIntersection (meshCenter, crP1, edge1ADT);
                bool p2e1 = checkEdgeIntersection (meshCenter, crP2, edge1ADT);

                bool p1e0 = checkEdgeIntersection (meshCenter, crP1, edge0ADT);
                bool p2e0 = checkEdgeIntersection (meshCenter, crP2, edge0ADT);*/
                
                /*int p1e1 = checkNumberOfEdgeIntersection (farPoint, crP1, edge1ADT);
                int p2e1 = checkNumberOfEdgeIntersection (farPoint, crP2, edge1ADT);
                
                int p1e0 = checkNumberOfEdgeIntersection (farPoint, crP1, edge0ADT);
                int p2e0 = checkNumberOfEdgeIntersection (farPoint, crP2, edge0ADT);*/
                
                
                bool crP1InsideDomain = rayCasting (farPoint, crP1, edgeADT);
                bool crP2InsideDomain = rayCasting (farPoint, crP2, edgeADT);
                
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
                    //cout << "intersects both in useNewPoints(...)" << endl;
                    cout << "both normal points are outside domain in useNewPoints(...)" << endl;
                    frontFirst.newPointChecked = true;
                    newPointSuccess = false;
                    points.pop_back(); // erase crP1
                    points.pop_back(); // erase crP2
                    //exit(-2);
                    return;
                }
                
                /*int p1e = checkNumberOfEdgeIntersection (farPoint, crP1, edgeADT);
                int p2e = checkNumberOfEdgeIntersection (farPoint, crP2, edgeADT);
                
                if (p1e % 2 != 0) // odd (inside)
                {
                    crP = crP1;
                }
                else if (p2e % 2 != 0)
                {
                    crP = crP2;
                }
                else
                {
                    //cout << "intersects both in useNewPoints(...)" << endl;
                    cout << "both normal points are outside domain in useNewPoints(...)" << endl;
                    frontFirst.newPointChecked = true;
                    newPointSuccess = false;
                    points.pop_back(); // erase crP1
                    points.pop_back(); // erase crP2
                    //exit(-2);
                    return;
                }*/

                /*if (frontEdge.belonging == 0)
                {
                    if (p1e0 % 2 != 0) // odd (inside)
                    {
                        crP = crP1;
                    }
                    else if (p2e0 % 2 != 0)
                    {
                        crP = crP2;
                    }
                    else
                    {
                        //cout << "intersects both in useNewPoints(...)" << endl;
                        cout << "both normal points are outside domain of mesh0 in useNewPoints(...)" << endl;
                        frontFirst.newPointChecked = true;
                        newPointSuccess = false;
                        points.pop_back(); // erase crP1
                        points.pop_back(); // erase crP2
                        //exit(-2);
                        return;
                    }
                }
                else if (frontEdge.belonging == 1)
                {
                    if (p1e1 % 2 != 0)
                    {
                        crP = crP1;
                    }
                    else if (p2e1 % 2 != 0)
                    {
                        crP = crP2;
                    }
                    else
                    {
                        //cout << "not intersects both in useNewPoints(...)" << endl;
                        cout << "both normal points are outside domain of mesh1 in useNewPoints(...)" << endl;
                        frontFirst.newPointChecked = true;
                        newPointSuccess = false;
                        points.pop_back(); // erase crP1
                        points.pop_back(); // erase crP2
                        return;
                        //exit(-2);
                    }
                }
                else
                {
                    cout << "!!! Error: frontEdge.belonging != 0 or frontEdge.belonging != 1 in useNewPoints(...)" << endl;
                    exit(-2);
                }*/
            }
            
            points.pop_back(); // erase crP1
            points.pop_back(); // erase crP2
            
            //cout << "df0 = " << points.size() << endl;
            //cin.ignore();
            
            bool doesPointExist = pointExists (crP.dim, points);
    
            CVector meshDis;
            meshDis[0] = 0.5 * pdis;
            meshDis[1] = 0.5 * pdis;
            meshDis[2] = 0.;
    
            bool pointNearby = pointsNearby (crP.dim-meshDis, crP.dim+meshDis, pointADT);
            
            if (!pointNearby)
            {
                if (!doesPointExist)
                {
                    bool intersection = checkEdgeIntersection (crP, points[it0], edgeADT);

                    if (!intersection)
                    {
                        intersection = checkEdgeIntersection (crP, points[it1], edgeADT);

                        if (!intersection)
                        {
                            crP.belonging = newGridId;
                            crP.newlyCreated = true;
                            
                            points.push_back (crP);
                            /*vec = pointADT.createADTPoint (crP.dim, crP.dim);
                            vec.idx = points.size() - 1;
                            pointADT.insert (vec, pointADT.root, tempBool);*/
                            
                            int iCrP = points.size() - 1;
                            
                            Edge tmpEdge1 = createEdge (iCrP, it0, newGridId, true);
                            Edge tmpEdge2 = createEdge (iCrP, it1, newGridId, true);
                            
                            edges.push_back (tmpEdge1);
                            edges.push_back (tmpEdge2);
                            
                            int iTmpEdge1 = edges.size() - 2;
                            int iTmpEdge2 = edges.size() - 1;
                                
                            tmpTriangle = createTriangle (iFrontEdge, iTmpEdge1, iTmpEdge2, edges, points);
                            bool triIntersects = triangleIntersect (tmpTriangle, triangleADT, points);

                            if (!triIntersects)
                            {
                                vec = pointADT.createADTPoint (crP.dim, crP.dim);
                                vec.idx = points.size() - 1;
                                pointADT.insert (vec, pointADT.root, tempBool);

                                vec = edgeADT.createADTPoint (crP, points[it0]);
                                vec.idx = iTmpEdge1;
                                edgeADT.insert (vec, edgeADT.root, tempBool);

                                vec = edgeADT.createADTPoint (crP, points[it1]);
                                vec.idx = iTmpEdge2;
                                edgeADT.insert (vec, edgeADT.root, tempBool);

                                triangles.push_back(tmpTriangle);

                                vec = triangleADT.createADTPoint (tmpTriangle, points);
                                vec.idx = triangles.size() - 1;
                                triangleADT.insert (vec, triangleADT.root, tempBool);

                                eraseFromFrontList (frontList);
                                addToFrontList (iTmpEdge2, frontList);
                                addToFrontList (iTmpEdge1, frontList);
                                sortFrontList (frontList, points, edges);

                                newPointSuccess = true;
                                return;
                            }
                            else
                            {
                                points.pop_back(); // erase crP
                                edges.pop_back(); // erase tmpEdge1
                                edges.pop_back(); // erase tmpEdge2
                            }
                        }
                    }
                }
            }
        }

        frontFirst.newPointChecked = true;
        newPointSuccess = false;
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
    
    bool rayCasting (const Point& farPoint, const Point& p, EdgeADT& edgeADT)
    {
        // determines whether point, p is inside domain of edgeADT or not
        // returns true if inside
        // returns false if outside
        
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
}


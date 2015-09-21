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
            srchCandPts (frontFirst, edges, points, pointADT, candPts, (2.*rho), edgeADT, edge01ADT);
            
            if (candPts.size() == 0)
            {
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
                construct (iCPX, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, it0,
                           it1, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points, edgeCenters, edgeCenterADT, circleADT);
            }
            else // if existing point doesn't work
            {
                // get a point normal to front
                Point crP = getNewPt (t0, t1, aveTriArea, edge01ADT);
                
                // no need to check circumradius for uniform grid
                // but will be needed in future
                
                // check two forming edges intersections
                bool newPointPass = checkTwoFormingEdges (crP, t0, t1, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, edges, edgeADT, points);
                
                if (newPointPass)
                {
                    // check if new point intersects circumcircles
                    if (!ptCCInter (crP, circleADT))
                    {
                        // add new point to points list
                        addToPointList (crP, points, pointADT);
                        int iCrP = points.size() - 1;
                        
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

                        outputTrianglesVTK (points, triangles, "../out", "tri.vtk");

                        exit(-2);
                    }
                }
                else
                {
                    cout << "none of the ways work in AFT::advanceFront(...)" << endl;
                
                    cout << "it0 = " << it0 << endl;
                    cout << "it1 = " << it1 << endl;
                    cout << "iCPX = " << iCPX << endl;

                    outputTrianglesVTK (points, triangles, "../out", "tri.vtk");

                    exit(-2);
                }
            }
            
            /*elAB = eligible (candPts.back(), points[candPts.back()], false, it0, it1, aveTriArea, scoreAB, A_CPX_exists, B_CPX_exists, iA_CPX,
            //                iB_CPX, frontList, edges, edgeADT, edge01ADT, triangleADT, points, pointADT, edgeCenterADT, circleADT);

            A_CPX_exists_AB[0] = A_CPX_exists;
            B_CPX_exists_AB[0] = B_CPX_exists;
            iA_CPX_AB[0] = iA_CPX;
            iB_CPX_AB[0] = iB_CPX;
            scores[0] = scoreAB;
            //ids[0] = frontFirst.cloPts.front();
            ids[0] = candPts.back();
            
            //
            double aveEdgeSize = sqrt ( (4./sqrt(3.) * aveTriArea) );
            getTwoNormalPoints (it0, it1, points, cnp1, cnp2, pdis, aveEdgeSize);
            points.push_back (cnp1); iCnp1 = points.size() - 1;
            points.push_back (cnp2); iCnp2 = points.size() - 1;
            
            
            
            elNP1 = eligible (iCnp1, points[iCnp1], true, it0, it1, aveTriArea, scoreNP1, A_CPX_exists, B_CPX_exists, iA_CPX,
                              iB_CPX, frontList, edges, edgeADT, edge01ADT, triangleADT, points, pointADT, edgeCenterADT, circleADT);
            
            A_CPX_exists_NP1 = A_CPX_exists;
            B_CPX_exists_NP1 = B_CPX_exists;
            iA_CPX_NP1 = iA_CPX;
            iB_CPX_NP1 = iB_CPX;
            
            elNP2 = eligible (iCnp2, points[iCnp2], true, it0, it1, aveTriArea, scoreNP2, A_CPX_exists, B_CPX_exists, iA_CPX,
                              iB_CPX, frontList, edges, edgeADT, edge01ADT, triangleADT, points, pointADT, edgeCenterADT, circleADT);
            
            A_CPX_exists_NP2 = A_CPX_exists;
            B_CPX_exists_NP2 = B_CPX_exists;
            iA_CPX_NP2 = iA_CPX;
            iB_CPX_NP2 = iB_CPX;
            
            cnp1 = points[iCnp1];
            cnp2 = points[iCnp2];
            
            points.pop_back();
            points.pop_back();
            //
            
            score = BIG_POS_NUM;
            if (elAB) { score = min (score, scoreAB); }
            if (elNP1) { score = min (score, scoreNP1); }
            if (elNP2) { score = min (score, scoreNP2); }
            
            if (elAB && elNP1)
            {
                if ( (scoreAB-scoreNP1) < tolPreferExistingPoint )
                {
                    score = scoreAB;
                }
            }
            
            if (elAB && elNP2)
            {
                if ( (scoreAB-scoreNP2) < tolPreferExistingPoint )
                {
                    score = scoreAB;
                }
            }
            
            
            
            if (score == BIG_POS_NUM)
            {
                //cout << "frontFirst.cloPts.size()1 = " << frontFirst.cloPts.size() << endl;
                //cout << "frontFirst.cloPts[0]()1 = " << frontFirst.cloPts[0] << endl;
                //frontFirst.cloPts.pop_front();
                int initSize = candPts.size();
                candPts.pop_back();
                
                //while (frontFirst.cloPts.size() != 0)
                while (candPts.size() != 0)
                {
                    //cout << "frontFirst.cloPts.size() = " << frontFirst.cloPts.size() << endl;
                    //cout << "frontFirst.cloPts[0] = " << frontFirst.cloPts[0] << endl;
                    
                    elAB = eligible (candPts.back(), points[candPts.back()], false, it0, it1, aveTriArea, scoreAB, A_CPX_exists, B_CPX_exists, iA_CPX,
                                     iB_CPX, frontList, edges, edgeADT, edge01ADT, triangleADT, points, pointADT, edgeCenterADT, circleADT);
                                        
                    //int i = frontFirst.cloPtsMaxSize - frontFirst.cloPts.size();
                    int i = initSize - candPts.size();
                    A_CPX_exists_AB[i] = A_CPX_exists;
                    B_CPX_exists_AB[i] = B_CPX_exists;
                    iA_CPX_AB[i] = iA_CPX;
                    iB_CPX_AB[i] = iB_CPX;
                    scores[i] = scoreAB;
                    //ids[i] = frontFirst.cloPts.front();
                    ids[i] = candPts.back();
                    
                    
                    
                    if (elAB)
                    {
                        //iChosenPoint = frontFirst.cloPts.front();
                        iChosenPoint = candPts.back();
                        isNewPoint = false;
                        break;
                    }
                    else
                    {
                        //frontFirst.cloPts.pop_front();
                        candPts.pop_back();
                    }
                    
                    //cout << "cloPts.size = " << frontFirst.cloPts.size() << endl;
                }
                
                if (candPts.size() == 0)
                {
                    double tmpScore = BIG_POS_NUM;                    
                    for (int i=0; i<initSize; ++i)
                    {
                        if (scores[i] < tmpScore)
                        {
                            iChosenPoint = ids[i];
                            isNewPoint = false;
                            elAB = true;
                            A_CPX_exists = A_CPX_exists_AB[i];
                            B_CPX_exists = B_CPX_exists_AB[i];
                            iA_CPX = iA_CPX_AB[i];
                            iB_CPX = iB_CPX_AB[i];
                            tmpScore = scores[i];
                        }
                    }
                    
                    if (!elAB)
                    {
                        cout << "none of the ways work in AFT::advanceFront(...)" << endl;
                    
                        cout << "score = " << score << endl;
                        cout << "scoreAB = " << scoreAB << endl;
                        cout << "scoreNP1 = " << scoreNP1 << endl;
                        cout << "scoreNP2 = " << scoreNP2 << endl;

                        cout << "t0 ID = " << it0 << endl;
                        cout << "t1 ID = " << it1 << endl;

                        cout << "t0Belo = " << points[it0].belonging << endl;
                        cout << "t1Belo = " << points[it1].belonging << endl;

                        cout << "t0.dim[0] = " << points[it0].dim[0] << endl;
                        cout << "t0.dim[1] = " << points[it0].dim[1] << endl;

                        cout << "t1.dim[0] = " << points[it1].dim[0] << endl;
                        cout << "t1.dim[1] = " << points[it1].dim[1] << endl;

                        cout << "triangles.size() = " << triangles.size() << endl;

                        cout << "edgeID = " << frontFirst.edge << endl;
                        cout << "edges.size() = " << edges.size() << endl;

                        cout << "newlyCreated = " << frontEdge.newlyCreated << endl;

                        //outputTriangles (points, triangles);
                        outputTrianglesVTK (points, triangles, "../out", "tri.vtk");

                        exit(-2);
                    }
                }
            }
            else if (score == scoreAB)
            {
                //iChosenPoint = frontFirst.cloPts.front();
                iChosenPoint = candPts.back();
                isNewPoint = false;
                A_CPX_exists = A_CPX_exists_AB[0];
                B_CPX_exists = B_CPX_exists_AB[0];
                iA_CPX = iA_CPX_AB[0];
                iB_CPX = iB_CPX_AB[0];
            }
            else if (score == scoreNP1)
            {
                addToPointList (cnp1, points, pointADT);
                iCnp1 = points.size() - 1;
                iChosenPoint = iCnp1;
                isNewPoint = true;
                A_CPX_exists = A_CPX_exists_NP1;
                B_CPX_exists = B_CPX_exists_NP1;
                iA_CPX = iA_CPX_NP1;
                iB_CPX = iB_CPX_NP1;
            }
            else if (score == scoreNP2)
            {
                addToPointList (cnp2, points, pointADT);
                iCnp2 = points.size() - 1;
                iChosenPoint = iCnp2;
                isNewPoint = true;
                A_CPX_exists = A_CPX_exists_NP2;
                B_CPX_exists = B_CPX_exists_NP2;
                iA_CPX = iA_CPX_NP2;
                iB_CPX = iB_CPX_NP2;
            }
            
            if (elAB || elNP1 || elNP2)
            {
                
                
                construct (iChosenPoint, isNewPoint, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, it0,
                           it1, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points, edgeCenters, edgeCenterADT, circleADT);
                
                
                cout << "elAB = " << elAB << endl;
                cout << "elNP1 = " << elNP1 << endl;
                cout << "elNP2 = " << elNP2 << endl;
                
                cout << "scoreAB = " << scoreAB << endl;
                cout << "scoreNP1 = " << scoreNP1 << endl;
                cout << "scoreNP2 = " << scoreNP2 << endl;
                
                exit(-2);
                
                
            }*/
            
            cout << "frontListSize = " << frontList.size() << endl;
        }
    }
}

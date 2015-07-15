#include "AFT.h"

namespace AFT
{
    FrontMember::FrontMember()
    {
        newPointChecked = false;
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
                       PointADT& pointADT, EdgeADT& edgeADT, EdgeADT& edge01ADT, int newGridId)
    {
        //This is specigic to 2 grid arrangement !!!!!!!!!!!!!!!!!!!!!
        
        while (!frontList.empty())
        {
            double scoreA = BIG_POS_NUM;
            double scoreB = BIG_POS_NUM;
            double scoreNP1 = BIG_POS_NUM;
            double scoreNP2 = BIG_POS_NUM;
            double score = BIG_POS_NUM;
            //bool CP0Checked = false;
            //bool CP1Checked = false;
            bool pointFound;
            //bool t0_t1_CP0_intersects;
            //bool t0_t1_CP1_intersects;
            bool A_CPX_exists;
            bool B_CPX_exists;
            int iA_CPX;
            int iB_CPX;
            bool elA;
            bool elB;
            bool elNP1;
            bool elNP2;
            Point cnp1;
            Point cnp2;
            int iCnp1;
            int iCnp2;
            int iChosenPoint;
            bool isNewPoint;
            double pdis;
            
            FrontMember& frontFirst = frontList.front();
            Edge& frontEdge = edges [ frontFirst.edge ];
            int it0 = frontEdge.t[0];
            int it1 = frontEdge.t[1];
            //const Point& t0 = points[ it0 ];
            //const Point& t1 = points[ it1 ];
            
            frontFirst.CPfound = 0;
            
            int iCP0 = findClosestPoint (frontFirst, 0, points, edges, pointFound);            
            if (pointFound)
            {
                elA = eligible (iCP0, false, it0, it1, aveTriArea, scoreA, A_CPX_exists, B_CPX_exists, iA_CPX,
                                iB_CPX, frontList, edges, edgeADT, edge01ADT, triangleADT, points, pointADT);
            }
            
            int iCP1 = findClosestPoint (frontFirst, 1, points, edges, pointFound);
            /*if (frontFirst.CPfound == 0)
            {                
                cout << "no close points found to each terminal" << endl;
                cout << t0.dim[0] << endl;
                cout << t0.dim[1] << endl;
                cout << t0.dim[2] << endl;
                
                cout << t1.dim[0] << endl;
                cout << t1.dim[1] << endl;
                cout << t1.dim[2] << endl;
                
                cout << triangles.size() << endl;
                
                cout << "belonging = " << frontEdge.belonging << endl;
                cout << "newlyCreated = " << frontEdge.newlyCreated << endl;
                cout << "t0 ID = " << it0 << endl;
                cout << "t1 ID = " << it1 << endl;
                
                //outputTriangles (points, triangles);
                outputTrianglesVTK (points, triangles, "../out");
                exit(-2);
            }*/            
            if (pointFound)
            {
                elB = eligible (iCP1, false, it0, it1, aveTriArea, scoreB, A_CPX_exists, B_CPX_exists, iA_CPX,
                                iB_CPX, frontList, edges, edgeADT, edge01ADT, triangleADT, points, pointADT);
            }              
            
            //
            getTwoNormalPoints (it0, it1, points, cnp1, cnp2, pdis);
            points.push_back (cnp1); iCnp1 = points.size() - 1;
            points.push_back (cnp2); iCnp2 = points.size() - 1;
            
            elNP1 = eligible (iCnp1, true, it0, it1, aveTriArea, scoreNP1, A_CPX_exists, B_CPX_exists, iA_CPX,
                              iB_CPX, frontList, edges, edgeADT, edge01ADT, triangleADT, points, pointADT);
            elNP2 = eligible (iCnp2, true, it0, it1, aveTriArea, scoreNP2, A_CPX_exists, B_CPX_exists, iA_CPX,
                              iB_CPX, frontList, edges, edgeADT, edge01ADT, triangleADT, points, pointADT);
            points.pop_back();
            points.pop_back();
            //
            
            if (elA) { score = min (score, scoreA); }
            if (elB) { score = min (score, scoreB); }
            if (elNP1) { score = min (score, scoreNP1); }
            if (elNP2) { score = min (score, scoreNP2); }
            
            if (score == BIG_POS_NUM)
            {
                cout << "none of the ways work in AFT::advanceFront(...)" << endl;
                cout << "score = " << score << endl;
                cout << "scoreA = " << scoreA << endl;
                cout << "scoreB = " << scoreB << endl;
                cout << "scoreNP1 = " << scoreNP1 << endl;
                cout << "scoreNP2 = " << scoreNP2 << endl;
                
                cout << "t0 ID = " << it0 << endl;
                cout << "t1 ID = " << it1 << endl;
                
                cout << "t0.dim[0] = " << points[it0].dim[0] << endl;
                cout << "t0.dim[1] = " << points[it0].dim[1] << endl;
                
                cout << "t1.dim[0] = " << points[it1].dim[0] << endl;
                cout << "t1.dim[1] = " << points[it1].dim[1] << endl;
                
                cout << "triangles.size() = " << triangles.size() << endl;
                
                cout << "belonging = " << frontEdge.belonging << endl;
                cout << "newlyCreated = " << frontEdge.newlyCreated << endl;
                
                //outputTriangles (points, triangles);
                outputTrianglesVTK (points, triangles, "../out");
                
                exit(-2);
            }
            else if (score == scoreA)
            {
                iChosenPoint = iCP0;
                isNewPoint = false;
            }
            else if (score == scoreB)
            {
                iChosenPoint = iCP1;
                isNewPoint = false;
            }
            else if (score == scoreNP1)
            {
                iChosenPoint = iCnp1;
                isNewPoint = true;
                points.push_back (cnp1);
                iCnp1 = points.size() - 1;
            }
            else if (score == scoreNP2)
            {
                iChosenPoint = iCnp2;
                isNewPoint = true;
                points.push_back (cnp2);
                iCnp2 = points.size() - 1;
            }
            
            construct (iChosenPoint, isNewPoint, A_CPX_exists, B_CPX_exists, iA_CPX, iB_CPX, it0,
                       it1, frontList, edges, triangles, triangleADT, edgeADT, newGridId, points);
                                    
            /*bool useNewPoint = false;
            //createNewPoint (useNewPoint, frontList, aveMeshSize, points, edges, triangles, triangleADT, pointADT, edgeADT, edge0ADT, edge1ADT, meshCenter, newGridId);            
                        
            if (!useNewPoint)
            {
                if (!CP0Checked || !CP1Checked)
                {
                    if (!CP0Checked && !CP1Checked)
                    {
                        if (disB < disA)
                        {
                            F1 (iCP1, iCP0, CP0Checked, CP1Checked, 1, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                        }
                        else
                        {
                            F1 (iCP0, iCP1, CP1Checked, CP0Checked, 0, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                        }
                    }
                    else if (!CP0Checked)
                    {
                        F1 (iCP0, iCP1, CP1Checked, CP0Checked, 0, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                    }
                    else
                    {
                        F1 (iCP1, iCP0, CP0Checked, CP1Checked, 1, frontList, edges, triangles, edgeADT, triangleADT, newGridId, points);
                    }
                }
            }*/
            
            cout << "frontListSize = " << frontList.size() << endl;
        }
    }
}
#include "AFT.h"

namespace AFT
{
    bool TriangleADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
    {
        bool insideCube = true;

        if (node->level != 0)
        {
            for (unsigned int d=0; d<ADT_DIM; ++d)
            {
                if (!(node->p.dim[d*2] <= targetPoint.dim[d*2+1]) || !(node->p.dim[d*2+1] >= targetPoint.dim[d*2]))
                {
                    insideCube = false;
                    break;
                }
            }
        }
        else
        {
            insideCube = false;
        }
        
        /*if ( fabs( targetPoint.dim[0] - 12.6295 ) < 0.01 )
        {
            if ( fabs (targetPoint.dim[1] - 16.106) < 0.01 )
            {
                cout << "found2" << endl;

                cout << node->p.dim[0] << endl;
                cout << node->p.dim[1] << endl;
                cout << node->p.dim[2] << endl;
                cout << node->p.dim[3] << endl;
                cout << node->p.dim[4] << endl;
                cout << node->p.dim[5] << endl;
                
                cout << targetPoint.dim[0] << endl;
                cout << targetPoint.dim[1] << endl;
                cout << targetPoint.dim[2] << endl;
                cout << targetPoint.dim[3] << endl;
                cout << targetPoint.dim[4] << endl;
                cout << targetPoint.dim[5] << endl;
                
                cin.ignore();
            }
        }*/

        return insideCube;
    }
    
    bool TriangleADT::compareFunction (const Node* node, const ADTPoint& targetPoint)
    {
        /*if ( fabs( targetPoint.dim[0] - 12.6295 ) < 0.01 )
        {
            if ( fabs (targetPoint.dim[1] - 16.106) < 0.01 )
            {
                cout << "found2" << endl;

                cout << node->p.dim[0] << endl;
                cout << node->p.dim[1] << endl;
                cout << node->p.dim[2] << endl;
                cout << node->p.dim[3] << endl;
                cout << node->p.dim[4] << endl;
                cout << node->p.dim[5] << endl;
                
                cout << targetPoint.dim[0] << endl;
                cout << targetPoint.dim[1] << endl;
                cout << targetPoint.dim[2] << endl;
                cout << targetPoint.dim[3] << endl;
                cout << targetPoint.dim[4] << endl;
                cout << targetPoint.dim[5] << endl;
                
                cin.ignore();
            }
        }*/
        
        /*if ( fabs( node->p.dim[0] - 12.6295 ) < 0.01 && fabs( targetPoint.dim[0] - 12.6295 ) < 0.01 )
        {
            if ( fabs (node->p.dim[1] - 16.0156) < 0.01 && fabs (targetPoint.dim[1] - 16.106) < 0.01 )
            {
                cout << "found2" << endl;

                cout << node->p.dim[0] << endl;
                cout << node->p.dim[1] << endl;
                cout << node->p.dim[2] << endl;
                cout << node->p.dim[3] << endl;
                cout << node->p.dim[4] << endl;
                cout << node->p.dim[5] << endl;
                
                cout << targetPoint.dim[0] << endl;
                cout << targetPoint.dim[1] << endl;
                cout << targetPoint.dim[2] << endl;
                cout << targetPoint.dim[3] << endl;
                cout << targetPoint.dim[4] << endl;
                cout << targetPoint.dim[5] << endl;
                
                cin.ignore();
            }
        }*/
        
        /* If two triangles are identical they would not be intersecting because whenever count=reqCount=2, the function will
         * return true without waiting for count=3.
         *
         * For the intersection of two triangles, an edge of one of the triangles should intersect two edges
         * of the other triangle so when count=reqCount=2 then function will return true.
         * 
         * If one triangle is inside the other one...
         */
        
        bool inside_1;
        bool inside_2;

        if (node->level != 0)
        {
            const Point& p1 = targetPoint.vertices[0];
            const Point& p2 = targetPoint.vertices[1];
            const Point& p3 = targetPoint.vertices[2];

            const Point& k1 = node->p.vertices[0];
            const Point& k2 = node->p.vertices[1];
            const Point& k3 = node->p.vertices[2];
            
            unsigned int count = 0;
            unsigned int reqCount = 2;

            bool dummyEM;
            
            if ( doIntersect (k1.dim, k2.dim, p1.dim, p2.dim, dummyEM) ) {++count;} if (count==reqCount) return true;
            if ( doIntersect (k1.dim, k2.dim, p1.dim, p3.dim, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k1.dim, k2.dim, p2.dim, p3.dim, dummyEM) ) {++count;} if (count==reqCount) return true;            

            if ( doIntersect (k1.dim, k3.dim, p1.dim, p2.dim, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k1.dim, k3.dim, p1.dim, p3.dim, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k1.dim, k3.dim, p2.dim, p3.dim, dummyEM) ) {++count;} if (count==reqCount) return true;            

            if ( doIntersect (k2.dim, k3.dim, p1.dim, p2.dim, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k2.dim, k3.dim, p1.dim, p3.dim, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k2.dim, k3.dim, p2.dim, p3.dim, dummyEM) ) {++count;} if (count==reqCount) return true;            

            // if one triangle is inside another one
            if (count <= 1) // count was equal to 1 before.
            {
                inside_1 = localCmpFunc (node->p, targetPoint);
                inside_2 = localCmpFunc (targetPoint, node->p);

                if (inside_1 || inside_2)
                {
                    return true;
                }
            }
        }
        else
        {
            return false;
        }

        return false;
    }
    
    bool TriangleADT::localCmpFunc (const ADTPoint& node, const ADTPoint& targetPoint)
    {
        bool inside;
        Point one;
        Point two;
        Point thr;
        Point cp1;
        Point cp2;        

        const Point&  p1 = targetPoint.vertices[0];
        const Point&  p2 = targetPoint.vertices[1];
        const Point&  p3 = targetPoint.vertices[2];

        const Point&  k1 = node.vertices[0];
        const Point&  k2 = node.vertices[1];
        const Point&  k3 = node.vertices[2];        

        for (unsigned int i=0; i<targetPoint.vertices.size(); ++i)
        {
            inside = true;

            one.dim = k1.dim - k2.dim;
            two.dim = targetPoint.vertices[i].dim - k2.dim;
            thr.dim = k3.dim - k2.dim;

            cp1.dim = crossP (one.dim, two.dim);
            cp2.dim = crossP (one.dim, thr.dim);

            if ( dotP (cp1.dim, cp2.dim) <= 0. ) inside = false;

            one.dim = k1.dim - k3.dim;
            two.dim = targetPoint.vertices[i].dim - k3.dim;
            thr.dim = k2.dim - k3.dim;

            cp1.dim = crossP (one.dim, two.dim);
            cp2.dim = crossP (one.dim, thr.dim);

            if ( dotP (cp1.dim, cp2.dim) <= 0. ) inside = false;

            one.dim = k2.dim - k3.dim;
            two.dim = targetPoint.vertices[i].dim - k3.dim;
            thr.dim = k1.dim - k3.dim;

            cp1.dim = crossP (one.dim, two.dim);
            cp2.dim = crossP (one.dim, thr.dim);

            if ( dotP (cp1.dim, cp2.dim) <= 0. ) inside = false;

            if (inside == true) return true;
        }

        return inside;
    }
    
    void TriangleADT::build (const EdgeADT& edgeADT)
    {
        this->root = new TriangleADT::Node();
        this->root->level = 0;
        this->root->p.idx = -1;

        for (unsigned int d=0; d<ADT_VAR; ++d)
        {
            this->root->c[d] = edgeADT.root->c[d];
            this->root->d[d] = edgeADT.root->d[d];            
        }
        
        

        // note that no point assigned to the root
    }
    
    ADT::ADTPoint TriangleADT::createADTPoint (const Triangle& tri, const vector<Point>& points)
    {
        ADTPoint vec;
        
        const Point& p0 = points[ tri.p[0] ];
        const Point& p1 = points[ tri.p[1] ];
        const Point& p2 = points[ tri.p[2] ];

        for (unsigned int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = min( min(p0.dim[i], p1.dim[i]), p2.dim[i] );
            vec.dim[i*2+1] = max( max(p0.dim[i], p1.dim[i]), p2.dim[i] );
        }

        vec.vertices.push_back (p0);
        vec.vertices.push_back (p1);
        vec.vertices.push_back (p2);

        return vec;
    }
    
    Triangle createTriangle (int e1, int e2, int e3, const vector<Edge>& edges, const vector<Point>& points)
    {
        Triangle tri;
        
        int e1t0 = edges[e1].t[0];
        int e1t1 = edges[e1].t[1];
        
        int e2t0 = edges[e2].t[0];
        int e2t1 = edges[e2].t[1];
        
        int e3t0 = edges[e3].t[0];
        int e3t1 = edges[e3].t[1];
        
        tri.e.push_back ( e1 );
        tri.e.push_back ( e2 );
        tri.e.push_back ( e3 );
        
        tri.p.push_back ( e1t0 );
        tri.p.push_back ( e1t1 );

        if ( (e2t0 != tri.p[0]) && (e2t0 != tri.p[1]) )
        {
            tri.p.push_back ( e2t0 );
        }
        else
        {
            tri.p.push_back ( e2t1 );
        }

        return tri;
    }
    
    bool triangleIntersect (const Triangle& tri, TriangleADT& triangleADT, const vector<Point>& points)
    {
        int result = -1;
        bool intersection = false;

        ADT::ADTPoint vec = triangleADT.createADTPoint (tri, points);

        triangleADT.searchForNIntersections = false;        
        result = triangleADT.search (vec);

        if (result != -1)
        {
            intersection = true;
        }

        return intersection;
    }
    
    double areaTriangle (const Triangle& tri, const vector<Point>& points)
    {
        const CVector& p0 = points[ tri.p[0] ].dim;
        const CVector& p1 = points[ tri.p[1] ].dim;
        const CVector& p2 = points[ tri.p[2] ].dim;

        CVector e0 = p1 - p0;
        CVector e1 = p2 - p0;

        return ( mag ( crossP(e0,e1) ) / 2. );
    }
    
    double charTriangleLength (const Triangle& tri, const vector<Point>& points)
    {
        return ( sqrt( areaTriangle(tri, points) ) );
    }
    
    void outputTriangles (const vector<Point>& points, const vector<Triangle>& triangles)
    {
        int len, maxLen = 20;

        ofstream out;
        out.open("../out/triangles.dat");

        out << "VARIABLES = ";
        out << "\"X\", "; // 1
        out << "\"Y\", "; // 2
        out << "\"Z\", "; // 3
        out << "\"I\", "; // 4
        out << std::endl;

        out << "ZONE ";
        out << "ZONETYPE=FETRIANGLE ";
        out << "N=";    
        out << points.size();
        out << " E=";
        out << triangles.size();
        out << " DATAPACKING=BLOCK ";
        out << "VARLOCATION=([4]=CELLCENTERED)";
        out << std::endl;

        // XYZ
        for (int d=0; d<N_DIM; ++d)
        {
            len = 0;
            for (const Point& p: points)
            {
                ++len;
                if (len == maxLen)
                {
                    out << endl;
                    len = 0;
                }

                out << p.dim[d] << ", ";
            }
            out << endl;
        }

        // I
        len = 0;
        for (unsigned int t=0; t<triangles.size(); ++t)
        {
            ++len;
            if (len == maxLen)
            {
                out << endl;
                len = 0;
            }

            out << t << ", ";
        }
        out << endl;

        // connectivity
        for (const Triangle& t: triangles)
        {
            for (const int p: t.p)
            {
                out << setw(10) << (&points[p] - &points.front()) + 1;
            }
            
            out << endl;
        }

        out.close();
    }
    
    void outputTrianglesVTK (const vector<Point>& points, const vector<Triangle>& triangles, string dir, string fileName)
    {
        int cellListSize = 0;

        string temps = fileName;
        //string temps = "tri.vtk";
        string slash = "/";
        dir.append (slash);
        dir.append (temps);

        ofstream out;
        out.open (dir);

        out << "# vtk DataFile Version 3.0" << endl;
        out << "Triangles" << endl;
        out << "ASCII" << endl;
        out << "DATASET UNSTRUCTURED_GRID" << endl;
        out << "POINTS " << points.size() << " float" << endl;

        for (const Point& p: points)
        {
            out << p.dim[0];
            out << " ";
            out << p.dim[1];
            out << " ";
            out << p.dim[2];
            out << endl;
        }

        // get cell list size
        for (const Triangle& t: triangles)
        {
            cellListSize += (t.p.size() + 1);
        }

        out << endl;    
        out << "CELLS " << triangles.size() << " " << cellListSize << endl;
        
        for (const Triangle& t: triangles)
        {
            out << t.p.size();
            out << " ";
            
            for (const int p: t.p)
            {
                out << &points[p] - &points.front();
                out << " ";
            }
            
            out << endl;
        }
        
        out << "CELL_TYPES " << triangles.size() << endl;
        for (unsigned int t=0; t<triangles.size(); ++t)
        {
            out << static_cast<int>(vtkCellType_t::TRI);
            out << endl;
        }
        
        out << "CELL_DATA " << triangles.size() << endl;
        
        out << "SCALARS " << "I " << "int " << "1" << endl;
        out << "LOOKUP_TABLE default" << endl;    
        for (unsigned int t=0; t<triangles.size(); ++t)
        {
            out << t << endl;
        }

        out.close();
    }
    
    void findNeighbors (const vector<Edge>& edges, vector<Triangle>& triangles)
    {
        for (unsigned int t=0; t<triangles.size(); ++t)
        {
            triangles[t].nei.resize (3);
            
            for (unsigned int i=0; i<triangles[t].e.size(); ++i)
            {
                const Edge& e = edges[ triangles[t].e[i] ];
                
                if ( e.nei.size() == 2 )
                {
                    if (e.nei[0] == t)
                    {
                        triangles[t].nei[i] = e.nei[1];
                    }
                    else if (e.nei[1] == t)
                    {
                        triangles[t].nei[i] = e.nei[0];
                    }
                    else
                    {
                        cout << "e.nei[x] != t in findNeighbors(...)" << endl;
                        exit(-2);
                    }
                }
                else
                {
                    cout << "e.nei.size() != 2 in findNeighbors(...)" << endl;
                    exit(-2);
                }
            }
        }
    }
    
    void flip (vector<Triangle>& triangles, vector<Edge>& edges, const vector<Point>& points)
    {
        CVector center;
        double radius;
        int iNei;
        int iExclusivePoint;
        int iOwnExclusivePoint;
        int iNeiEdge;
        int iEdgeA;
        int iEdgeB;
        int iEdgeAN;
        int iEdgeBN;
        double disExclusivePointCenter;

        // note that rectangular coordiantes are considered not circular

        for (unsigned int t=0; t<triangles.size(); ++t)
        {
            Triangle& tri = triangles[t];
            
            const int ip0 = tri.p[0];
            const int ip1 = tri.p[1];
            const int ip2 = tri.p[2];
            const CVector& pd0 = points[ ip0 ].dim;
            const CVector& pd1 = points[ ip1 ].dim;
            const CVector& pd2 = points[ ip2 ].dim;
            
            circleCenter(pd0, pd1, pd2, center, radius);

            if ( true )
            //if ( isfinite(center[0]) && isfinite(center[1]) )
            {
                for (unsigned int n=0; n<tri.nei.size(); ++n)
                {
                    iNei = tri.nei[n];
                    Triangle& nei = triangles[ iNei ];
                    
                    if (iNei != -1)
                    {
                        iNeiEdge = tri.e[n];
                        Edge& neiEdge = edges[ iNeiEdge ];
                        int iNeiEdgeT0 = neiEdge.t[0];
                        int iNeiEdgeT1 = neiEdge.t[1];

                        // get exclusive point
                        for (const int p: nei.p)
                        {
                            if (p != iNeiEdgeT0 && p != iNeiEdgeT1)
                            {
                                iExclusivePoint = p;
                                break;
                            }
                        }                        
                        const Point& exclusivePoint = points[ iExclusivePoint ];

                        disExclusivePointCenter = sqrt(
                                pow(exclusivePoint.dim[0] - center[0],2) + pow(exclusivePoint.dim[1] - center[1],2)
                                );

                        if (disExclusivePointCenter <= radius)
                        {
                            // get own exclusive point
                            for (const int p: tri.p)
                            {
                                if (p != iNeiEdgeT0 && p != iNeiEdgeT1)
                                {
                                    iOwnExclusivePoint = p;                                    
                                    break;
                                }
                            }                            
                            const Point& ownExclusivePoint = points[iOwnExclusivePoint];

                            // get edgeA and edgeB
                            for (int e: tri.e)
                            {
                                if (e != iNeiEdge)
                                {
                                    if (edges[e].t[0] != iNeiEdgeT1 && edges[e].t[1] != iNeiEdgeT1)
                                    {
                                        iEdgeA = e;
                                    }
                                    else if (edges[e].t[0] != iNeiEdgeT0 && edges[e].t[1] != iNeiEdgeT0)
                                    {
                                        iEdgeB = e;
                                    }
                                }
                            }

                            // get edgeAN and edgeBN
                            for (int e: nei.e)
                            {
                                if (e != iNeiEdge)
                                {
                                    if (edges[e].t[0] != iNeiEdgeT1 && edges[e].t[1] != iNeiEdgeT1)
                                    {
                                        iEdgeAN = e;
                                    }
                                    else if (edges[e].t[0] != iNeiEdgeT0 && edges[e].t[1] != iNeiEdgeT0)
                                    {
                                        iEdgeBN = e;
                                    }
                                }
                            }
                            
                            Edge& edgeA  = edges[ iEdgeA ];
                            Edge& edgeB  = edges[ iEdgeB ];
                            Edge& edgeAN = edges[ iEdgeAN ];
                            Edge& edgeBN = edges[ iEdgeBN ];

                            for (unsigned int i=0; i<edgeB.nei.size(); ++i)
                            {
                                if (edgeB.nei[i] == t)
                                {
                                    int iTriB = edgeB.nei[1-i];
                                    
                                    if ( edgeB.nei[1-i] != -1 )
                                    {
                                        for (int& r: triangles[iTriB].nei)
                                        {
                                            if (r == t)
                                            {
                                                r = iNei;
                                                break;
                                            }
                                        }
                                    }

                                    for (int r=0; r<3; ++r)
                                    {
                                        if (nei.e[r] == iEdgeAN)
                                        {
                                            nei.e[r] = iEdgeB;
                                            nei.nei[r] = iTriB;
                                            break;
                                        }
                                    }

                                    edgeB.nei[i] = iNei;
                                    break;
                                }
                            }

                            for (unsigned int i=0; i<edgeAN.nei.size(); ++i)
                            {
                                if (edgeAN.nei[i] == iNei)
                                {
                                    int iTriAN = edgeAN.nei[1-i];
                                    
                                    if ( edgeAN.nei[1-i] != -1 )
                                    {
                                        Triangle& triAN = triangles[ iTriAN ];

                                        for (int& r: triAN.nei)
                                        {
                                            if (r == iNei)
                                            {
                                                r = t;
                                                break;
                                            }
                                        }
                                    }
                                    
                                    for (int r=0; r<3; ++r)
                                    {
                                        if (tri.e[r] == iEdgeB)
                                        {
                                            tri.e[r] = iEdgeAN;
                                            tri.nei[r] = iTriAN;
                                            break;
                                        }
                                    }

                                    edgeAN.nei[i] = t;
                                    break;
                                }
                            }
                            
                            tri.p[0] = iNeiEdgeT0;
                            tri.p[1] = iExclusivePoint;
                            tri.p[2] = iOwnExclusivePoint;
                            
                            nei.p[0] = iNeiEdgeT1;
                            nei.p[1] = iExclusivePoint;
                            nei.p[2] = iOwnExclusivePoint;

                            neiEdge.t[0] = iExclusivePoint;
                            neiEdge.t[1] = iOwnExclusivePoint;
                        }
                    }
                }
            }
            else
            {
                cout << "not finite in flip(...)" << endl;
                exit(-2);
            }
        }
    }
    
    void circleCenter(const CVector& A, const CVector& B, const CVector& C, CVector& cnt, double& radius)
    {
        // http://stackoverflow.com/questions/4103405/what-is-the-algorithm-for-finding-the-center-of-a-circle-from-three-points

        double xDelta_a, yDelta_a;
        double xDelta_b, yDelta_b;

        yDelta_a = B[1] - A[1];
        xDelta_a = B[0] - A[0];

        if (xDelta_a == 0.)
        {
            yDelta_a = C[1] - A[1];
            xDelta_a = C[0] - A[0];

            yDelta_b = C[1] - B[1];
            xDelta_b = C[0] - B[0];

            if (xDelta_b == 0.)
            {
                yDelta_b = B[1] - A[1];
                xDelta_b = B[0] - A[0];
            }
        }
        else
        {
            yDelta_b = C[1] - B[1];
            xDelta_b = C[0] - B[0];

            if (xDelta_b == 0.)
            {
                yDelta_b = C[1] - A[1];
                xDelta_b = C[0] - A[0];
            }
        }

        double aSlope = yDelta_a/xDelta_a;
        double bSlope = yDelta_b/xDelta_b;

        if (xDelta_a == 0. || xDelta_b == 0.)
        {
            cout << "yesso: " << endl;
            cin.ignore();
        }

        if (aSlope == bSlope)
        {
            cout << "norto: " << endl;
            cin.ignore();
        }

        cnt[0] = (aSlope*bSlope*(A[1] - C[1]) + bSlope*(A[0] + B[0]) - aSlope*(B[0]+C[0])) / (2.* (bSlope-aSlope));

        if (aSlope != 0.)
        {
            cnt[1] = -1. * (cnt[0] - (A[0] + B[0]) / 2.) / aSlope +  (A[1] + B[1]) / 2.;
        }
        else if (bSlope != 0.)
        {
            cnt[1] = -1. * (cnt[0] - (B[0] + C[0]) / 2.) / bSlope + (B[1] + C[1]) / 2.;
        }
        else
        {
            cnt[1] = 0.; // not needed
            cout << "WTF" << endl;
            cin.ignore();
        }

        double difX = A[0] - cnt[0];
        double difY = A[1] - cnt[1];

        radius = sqrt(difX * difX + difY * difY);
    }
    
    void addToTriangleList(vector<Triangle>& triangles, const Triangle& tmpTriangle,
            TriangleADT& triangleADT, const vector<Point>& points)
    {
        bool tmpBool;
        
        triangles.push_back (tmpTriangle);
        ADT::ADTPoint vec = triangleADT.createADTPoint (tmpTriangle, points);
        vec.idx = triangles.size() - 1;
        triangleADT.insert (vec, triangleADT.root, tmpBool);
    }
    
    CVector Triangle::centroid(const vector<Point>& points)
    {
        CVector cent;
        for (int i=0; i<N_DIM; ++i)
        {
            cent[i] = (points[p[0]].dim[i] + points[p[1]].dim[i] + points[p[2]].dim[i]) / 3.;
        }
        
        return cent;
    }
    
    double Triangle::qualityScore (const vector<Point>& points, double aveTriArea, bool verbose, bool& passed)
    {
        // low score is better
        
        #include "Triangle.h"

        passed = true;
        
        double a = mag( points[p[0]].dim - points[p[1]].dim );
        double b = mag( points[p[0]].dim - points[p[2]].dim );
        double c = mag( points[p[1]].dim - points[p[2]].dim );
        
        double s = 0.5 * (a + b + c);        
        double r = (a*b*c) / (4. * sqrt( s*(s-a)*(s-b)*(s-c) ) );
        double areaEqui = 3*sqrt(3.)/4. * pow(r, 2.);
        double areaCurrent = areaTriangle (*this, points);
        
        // skewness
        double skew = (1. - areaCurrent/areaEqui);
        
        // area
        double areaSmall = min (areaCurrent, aveTriArea);
        double areaLarge = max (areaCurrent, aveTriArea);        
        double devAveArea = areaLarge / areaSmall;
        //double devAveArea = fabs(areaCurrent - aveTriArea) / aveTriArea;
        
        // aspect ratio
        double lon = max (max(a, b), c);
        double sho = min (min(a, b), c);
        //double aR = lon/sho;
        double aR = (lon-sho) / sho;
        
        double rSkew = cSkew*skew;
        double rArea = cAA*devAveArea;
        double rAR = cAR*aR;
        double score = rSkew + rArea + rAR;
        
        if (skew >= maxSkew)
        {
            if (verbose)
            {
                cout << "skew > maxSkew " << endl;
                cout << "skew = " << skew << endl;
                cout << "maxSkew = " << maxSkew << endl;
                cout << "areaCurrent = " << areaCurrent << endl;
                cout << "areaEqui = " << areaEqui << endl;                
            }            
            
            passed = false;
            return score;
        }
        
        if (devAveArea > maxValArea)
        {
            if (verbose)
            {
                cout << "devAveArea > maxValArea" << endl;
                cout << "devAveArea = " << devAveArea << endl;
                cout << "maxValArea = " << maxValArea << endl;
                cout << "areaCurrent = " << areaCurrent << endl;
                cout << "aveTriArea = " << aveTriArea << endl;                
            }
            
            passed = false;
            return score;
        }
        
        if (aR > maxValAR)
        {
            if (verbose)
            {
                cout << "aR > maxValAR" << endl;
                cout << "aR = " << aR << endl;
                cout << "maxValAR = " << maxValAR << endl;
                cout << "lon = " << lon << endl;
                cout << "sho = " << sho << endl;
            }
            
            passed = false;
            return score;
        }
        
        if (verbose)
        {
            cout << "skew = " << skew << endl;
            cout << "devAveArea = " << devAveArea << endl;
            cout << "aR = " << aR << endl;
            cout << "rSkew = " << rSkew << endl;
            cout << "rArea = " << rArea << endl;
            cout << "rAR = " << rAR << endl;
            cout << "score = " << score << endl;
        }
        
        return score;
    }
}


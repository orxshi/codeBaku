#include "AFT.h"

namespace AFT
{
    Triangle::Triangle ()
    {
        alive = true;
    }

    void CircleADT::build (const EdgeADT& edgeADT)
    {
        this->root = new CircleADT::Node();
        this->root->level = 0;
        //this->root->p->idx = -1;
        idsInTree.push_back (-1);
        addrsInTree.push_back (root);

        for (unsigned int d=0; d<ADT_VAR; ++d)
        {
            this->root->c[d] = edgeADT.root->c[d];
            this->root->d[d] = edgeADT.root->d[d];            
        }

        // note that no point assigned to the root
    }
    
    bool CircleADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
    {
        bool insideCube = true;

        if (node->level != 0)
        {
            for (unsigned int d=0; d<ADT_DIM; ++d)
            {
                if (!(node->p->dim[d*2] <= targetPoint.dim[d*2+1]) || !(node->p->dim[d*2+1] >= targetPoint.dim[d*2]))
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

        return insideCube;
    }
    
    bool CircleADT::compareFunction (const Node* node, const ADTPoint& targetPoint)
    {
        // check whether one of the corners of rectangle intersects circle
        
        /*
         * 0: xmin
         * 1: xmax
         * 2: ymin
         * 3: ymax
         * 4: zmin
         * 5: zmax
         */
        
        // vectors of 4 corners of rectangle
        CVector corner[4];
        
        // lower left
        corner[0][0] = targetPoint.dim[0]; // xmin
        corner[0][1] = targetPoint.dim[2]; // ymin
        
        // upper left
        corner[1][0] = targetPoint.dim[0]; // xmin
        corner[1][1] = targetPoint.dim[3]; // ymax
        
        // upper right
        corner[2][0] = targetPoint.dim[1]; // xmax
        corner[2][1] = targetPoint.dim[3]; // ymax
        
        // lower right
        corner[3][0] = targetPoint.dim[1]; // xmax
        corner[3][1] = targetPoint.dim[2]; // ymin
        
        // z-components are zero due to 2D
        corner[0][2] = 0.;
        corner[1][2] = 0.;
        corner[2][2] = 0.;
        corner[3][2] = 0.;
        
        // center of circle
        CVector center = node->p->vertices[0];
        
        // radius of circle
        double r = fabs( node->p->dim[0] - node->p->vertices[0][0] ); // |xminRec - xCen|
        
        for (int i=0; i<4; ++i)
        {
            CVector d = corner[i] - center;
            
            if (mag(d) < r)
            {
                // corner[i] is inside circle
                return true;
            }
        }
        
        return false;
    }
    
    ADT::ADTPoint CircleADT::createADTPoint (const Triangle& tri, const vector<Point>& points, int idTri)
    {
        const CVector& p0 = points[ tri.p[0] ].dim;
        const CVector& p1 = points[ tri.p[1] ].dim;
        const CVector& p2 = points[ tri.p[2] ].dim;
        
        CVector cnt;
        double radius;
        
        triPtsCircums (p0, p1, p2, cnt, radius);
        
        ADTPoint vec;
        
        vec.dim[0] = cnt[0] - radius; // xmin
        vec.dim[1] = cnt[0] + radius; // xmax
        
        vec.dim[2] = cnt[1] - radius; // ymin
        vec.dim[3] = cnt[1] + radius; // ymax
        
        vec.dim[4] = 0.; // zmin
        vec.dim[5] = 0.; // zmax

        Point p;
        p.dim = cnt;
        
        vec.idx = idTri;        
        
        vec.vertices.push_back (p.dim);

        return vec;
    }
    
    bool TriangleADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
    {
        bool insideCube = true;

        if (node->level != 0)
        {
            for (unsigned int d=0; d<ADT_DIM; ++d)
            {
                if (!(node->p->dim[d*2] <= targetPoint.dim[d*2+1]) || !(node->p->dim[d*2+1] >= targetPoint.dim[d*2]))
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

        return insideCube;
    }
    
    bool TriangleADT::compareFunction (const Node* node, const ADTPoint& targetPoint)
    {
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
            const CVector& p1 = targetPoint.vertices[0];
            const CVector& p2 = targetPoint.vertices[1];
            const CVector& p3 = targetPoint.vertices[2];

            const CVector& k1 = node->p->vertices[0];
            const CVector& k2 = node->p->vertices[1];
            const CVector& k3 = node->p->vertices[2];
            
            unsigned int count = 0;
            unsigned int reqCount = 2;

            bool dummyEM;
            
            if ( doIntersect (k1, k2, p1, p2, dummyEM) ) {++count;} if (count==reqCount) return true;
            if ( doIntersect (k1, k2, p1, p3, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k1, k2, p2, p3, dummyEM) ) {++count;} if (count==reqCount) return true;            

            if ( doIntersect (k1, k3, p1, p2, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k1, k3, p1, p3, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k1, k3, p2, p3, dummyEM) ) {++count;} if (count==reqCount) return true;            

            if ( doIntersect (k2, k3, p1, p2, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k2, k3, p1, p3, dummyEM) ) {++count;} if (count==reqCount) return true;            
            if ( doIntersect (k2, k3, p2, p3, dummyEM) ) {++count;} if (count==reqCount) return true;            

            // if one triangle is inside another one
            if (count <= 1) // count was equal to 1 before.
            {
                inside_1 = localCmpFunc (*(node->p), targetPoint);
                inside_2 = localCmpFunc (targetPoint, *(node->p));

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
        CVector one;
        CVector two;
        CVector thr;
        CVector cp1;
        CVector cp2;        

        const CVector&  p1 = targetPoint.vertices[0];
        const CVector&  p2 = targetPoint.vertices[1];
        const CVector&  p3 = targetPoint.vertices[2];

        const CVector&  k1 = node.vertices[0];
        const CVector&  k2 = node.vertices[1];
        const CVector&  k3 = node.vertices[2];        

        for (unsigned int i=0; i<targetPoint.vertices.size(); ++i)
        {
            inside = true;

            one = k1 - k2;
            two = targetPoint.vertices[i] - k2;
            thr = k3 - k2;

            cp1 = crossP (one, two);
            cp2 = crossP (one, thr);

            if ( dotP (cp1, cp2) <= 0. ) inside = false;

            one = k1 - k3;
            two = targetPoint.vertices[i] - k3;
            thr = k2 - k3;

            cp1 = crossP (one, two);
            cp2 = crossP (one, thr);

            if ( dotP (cp1, cp2) <= 0. ) inside = false;

            one = k2 - k3;
            two = targetPoint.vertices[i] - k3;
            thr = k1 - k3;

            cp1 = crossP (one, two);
            cp2 = crossP (one, thr);

            if ( dotP (cp1, cp2) <= 0. ) inside = false;

            if (inside == true) return true;
        }

        return inside;
    }
    
    void TriangleADT::build (const EdgeADT& edgeADT)
    {
        this->root = new TriangleADT::Node();
        this->root->level = 0;
        //this->root->p->idx = -1;
        
        idsInTree.push_back (-1);
        addrsInTree.push_back (root);

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
        
        const CVector& p0 = points[ tri.p[0] ].dim;
        const CVector& p1 = points[ tri.p[1] ].dim;
        const CVector& p2 = points[ tri.p[2] ].dim;

        for (unsigned int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = min( min(p0[i], p1[i]), p2[i] );
            vec.dim[i*2+1] = max( max(p0[i], p1[i]), p2[i] );
        }

        vec.vertices.push_back (p0);
        vec.vertices.push_back (p1);
        vec.vertices.push_back (p2);

        return vec;
    }
    
    ADT::ADTPoint TriangleADT::createADTPoint (const Point& p0, const Point& p1, const Point& p2)
    {
        ADTPoint vec;
        
        //const Point& p0 = points[ tri.p[0] ];
        //const Point& p1 = points[ tri.p[1] ];
        //const Point& p2 = points[ tri.p[2] ];

        for (unsigned int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = min( min(p0.dim[i], p1.dim[i]), p2.dim[i] );
            vec.dim[i*2+1] = max( max(p0.dim[i], p1.dim[i]), p2.dim[i] );
        }

        vec.vertices.push_back (p0.dim);
        vec.vertices.push_back (p1.dim);
        vec.vertices.push_back (p2.dim);

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
    
    bool triangleIntersect (const Point& p0, const Point& p1, const Point& p2, TriangleADT& triangleADT)
    {
        int result = -1;
        bool intersection = false;

        //ADT::ADTPoint vec = triangleADT.createADTPoint (tri, points);
        ADT::ADTPoint vec = triangleADT.createADTPoint (p0, p1, p2);

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
        /*for (unsigned int t=0; t<triangles.size(); ++t)
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
        }*/
        
        for (const Triangle& t: triangles)
        {
            if (t.nei.size() != 2)
            {
                cout << "triangle nei size must be 2 in AFT::findNeighbors(...)" << endl;
                exit(-2);
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

        /*cout << "triangles[171].p[0] = " << triangles[171].p[0] << endl;
        cout << "triangles[171].p[1] = " << triangles[171].p[1] << endl;
        cout << "triangles[171].p[2] = " << triangles[171].p[2] << endl;*/
        
        for (unsigned int t=0; t<triangles.size(); ++t)
            //for (unsigned int t=0; t<170; ++t)
        {
            Triangle& tri = triangles[t];
            
            const int ip0 = tri.p[0];
            const int ip1 = tri.p[1];
            const int ip2 = tri.p[2];            
            const CVector& pd0 = points[ ip0 ].dim;
            const CVector& pd1 = points[ ip1 ].dim;
            const CVector& pd2 = points[ ip2 ].dim;
            
            /*if (t == 3)
            {
                cout << "tri.p[0] = " << tri.p[0] << endl;
                cout << "tri.p[1] = " << tri.p[1] << endl;
                cout << "tri.p[2] = " << tri.p[2] << endl;
                
            }
            
            if (tri.p[0] == tri.p[1] || tri.p[0] == tri.p[2] || tri.p[1] == tri.p[2])
            {
                cout << "points are same in flip triangles 1" << endl;
                cout << "tri.p[0] = " << tri.p[0] << endl;
                cout << "tri.p[1] = " << tri.p[1] << endl;
                cout << "tri.p[2] = " << tri.p[2] << endl;
                cout << "t = " << t << endl;
                exit(-2);
            }*/
            
            
            
            triPtsCircums (pd0, pd1, pd2, center, radius);
            
            /*if (t == 3)
                        {
                            exit(-2);
                        }*/

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
                        
                        /*if (t == 3)
                        {
                            cout << "radius = " << radius << endl;
                        }*/
                        
                        
                        if (disExclusivePointCenter <= radius)
                        {
                            
                            /*if (t == 3)
                            {
                                cout << "radius = " << radius << endl;
                                cout << "tri.p[0] = " << tri.p[0] << endl;
                                cout << "tri.p[1] = " << tri.p[1] << endl;
                                cout << "tri.p[2] = " << tri.p[2] << endl;
                                cout << "pd0[0] = " << pd0[0] << endl;
                                cout << "pd0[1] = " << pd0[1] << endl;
                                cout << "pd1[0] = " << pd1[0] << endl;
                                cout << "pd1[1] = " << pd1[1] << endl;
                                cout << "pd2[0] = " << pd2[0] << endl;
                                cout << "pd2[1] = " << pd2[1] << endl;
                                cout << "t = " << t << endl;
                                exit(-2);
                            }*/
                            
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
                            
                            /*if (t == 3)
                            {
                                cout << "tri.p[0] = " << tri.p[0] << endl;
                                cout << "tri.p[1] = " << tri.p[1] << endl;
                                cout << "tri.p[2] = " << tri.p[2] << endl;
                                
                            }
                            
                            
                            
                            if (tri.p[0] == tri.p[1] || tri.p[0] == tri.p[2] || tri.p[1] == tri.p[2])
                            {
                                cout << "points are same in flip triangles" << endl;
                                cout << "tri.p[0] = " << tri.p[0] << endl;
                                cout << "tri.p[1] = " << tri.p[1] << endl;
                                cout << "tri.p[2] = " << tri.p[2] << endl;
                                cout << "t = " << t << endl;
                                exit(-2);
                            }*/
                            
                            nei.p[0] = iNeiEdgeT1;
                            nei.p[1] = iExclusivePoint;
                            nei.p[2] = iOwnExclusivePoint;
                            
                            /*if (t == 3)
                            {
                                cout << "nei.p[0] = " << nei.p[0] << endl;
                                cout << "nei.p[1] = " << nei.p[1] << endl;
                                cout << "nei.p[2] = " << nei.p[2] << endl;
                                cout << "iNei = " << iNei << endl;
                                exit(-2);
                            }
                            
                            if (nei.p[0] == nei.p[1] || nei.p[0] == nei.p[2] || nei.p[1] == nei.p[2])
                            {
                                cout << "points are same in flip triangles for nei" << endl;
                                cout << "nei.p[0] = " << nei.p[0] << endl;
                                cout << "nei.p[1] = " << nei.p[1] << endl;
                                cout << "nei.p[2] = " << nei.p[2] << endl;
                                
                                exit(-2);
                            }*/

                            neiEdge.t[0] = iExclusivePoint;
                            neiEdge.t[1] = iOwnExclusivePoint;
                            
                            /*if (triangles[171].p[0] != 178 || triangles[171].p[2] != 14)
                            {
                                cout << t << endl;
                                cout << triangles[171].p[0] << endl;
                                cout << triangles[171].p[2] << endl;
                                exit(-2);
                            }*/
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
    
    double triEdgeCircumradius (double a, double b, double c)
    {
        return ( (a * b * c) / sqrt( (a+b+c)*(b+c-a)*(c+a-b)*(a+b-c) ) );
    }
    
    void triPtsCircums (const CVector& A, const CVector& B, const CVector& C, CVector& cnt, double& radius)
    {
        // http://paulbourke.net/geometry/circlesphere/

        double xDelta_a, yDelta_a;
        double xDelta_b, yDelta_b;
        double x1, x2, x3;
        double y1, y2, y3;
        
        x1 = A[0];        
        y1 = A[1];
        
        x2 = B[0];        
        y2 = B[1];
        
        x3 = C[0];
        y3 = C[1];

        if (fabs(x2 - x1) < 1e-5)
        {
            // swap p2 and p3
            
            x2 = C[0];            
            y2 = C[1];
            
            x3 = B[0];            
            y3 = B[1];
        }
        else
        {
            if (fabs(x2 - x3) < 1e-5)
            {
                // swap p1 and p2
                
                x1 = B[0];            
                y1 = B[1];

                x2 = A[0];        
                y2 = A[1];
            }
        }
        
        xDelta_a = x2 - x1;
        yDelta_a = y2 - y1;
        
        xDelta_b = x3 - x2;
        yDelta_b = y3 - y2;

        double aSlope = yDelta_a / xDelta_a;
        double bSlope = yDelta_b / xDelta_b;

        if (xDelta_a == 0. || xDelta_b == 0.)
        {
            cout << "yesso" << endl;
            exit(-2);
            //cin.ignore();
        }

        if (aSlope == bSlope)
        {
            cout << "aSlope = bSlope" << endl;
            cout << "aSlope =" << aSlope << endl;
            cout << "yDelta_a =" << yDelta_a << endl;
            cout << "xDelta_a =" << xDelta_a << endl;
            cout << "yDelta_b =" << yDelta_b << endl;
            cout << "xDelta_b =" << xDelta_b << endl;
            cout << "A[0] =" << A[0] << endl;
            cout << "A[1] =" << A[1] << endl;
            cout << "B[0] =" << B[0] << endl;
            cout << "B[1] =" << B[1] << endl;
            cout << "C[0] =" << C[0] << endl;
            cout << "C[1] =" << C[1] << endl;
            exit(-2);
        }

        cnt[0] = (aSlope*bSlope*(y1 - y3) + bSlope*(x1 + x2) - aSlope*(x2 + x3)) / (2.* (bSlope-aSlope));

        if (aSlope != 0.)
        {
            cnt[1] = -1. * (cnt[0] - (x1 + x2) / 2.) / aSlope + (y1 + y2) / 2.;
        }
        else if (bSlope != 0.)
        {
            cnt[1] = -1. * (cnt[0] - (x2 + x3) / 2.) / bSlope + (y2 + y3) / 2.;
        }        

        double difX = x1 - cnt[0];
        double difY = y1 - cnt[1];

        radius = sqrt(difX * difX + difY * difY);
    }
    
    void addToTriangleList(vector<Triangle>& triangles, Triangle& tmpTriangle,
            TriangleADT& triangleADT, vector<Point>& points, CircleADT& circleADT, vector<Edge>& edges)
    {
        bool tmpBool;
        
        triangles.push_back (tmpTriangle);
        ADT::ADTPoint vec = triangleADT.createADTPoint (tmpTriangle, points);
        vec.idx = triangles.size() - 1;
        triangleADT.insert (vec, triangleADT.root, tmpBool);
        
        ADT::ADTPoint vec1 = circleADT.createADTPoint (tmpTriangle, points, vec.idx);        
        circleADT.insert (vec1, circleADT.root, tmpBool);
        
        points[tmpTriangle.p[0]].tri.push_back (triangles.size() - 1);
        points[tmpTriangle.p[1]].tri.push_back (triangles.size() - 1);
        points[tmpTriangle.p[2]].tri.push_back (triangles.size() - 1);
        
        edges[tmpTriangle.e[0]].tri.push_back (triangles.size() - 1);
        edges[tmpTriangle.e[1]].tri.push_back (triangles.size() - 1);
        edges[tmpTriangle.e[2]].tri.push_back (triangles.size() - 1);
        
        for (int e: tmpTriangle.e)
        {
            if (edges[e].nei.size() > 1)
            {
                cout << "edges[e].nei.size cannot be greater than 1 in AFT::addToTriangleList(...)" << endl;
                
                for (int i=0; i<edges[e].nei.size(); ++i)
                {
                    cout << "edges[e].nei = " << edges[e].nei[i] << endl;
                }
                
                for (int i=0; i<tmpTriangle.e.size(); ++i)
                {
                    cout << "tmpTriangle.e = " << tmpTriangle.e[i] << endl;
                }
                
                cout << "t = " << triangles.size() - 1 << endl;
                cout << "e = " << e << endl;
                cout << "t[0] = " << edges[e].t[0] << endl;
                cout << "t[1] = " << edges[e].t[1] << endl;
                
                cout << "p[0] = " << tmpTriangle.p[0] << endl;
                cout << "p[1] = " << tmpTriangle.p[1] << endl;
                cout << "p[2] = " << tmpTriangle.p[2] << endl;
                
                outputTrianglesVTK (points, triangles, "../out", "tri.vtk");
                
                
                exit(-2);
            }
            
            if (edges[e].nei.size() == 1)
            {
                tmpTriangle.nei.push_back (edges[e].nei[0]); // you need do this before the following statement
            }
            
            edges[e].nei.push_back (triangles.size() - 1);
        }
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
    
    CVector cntTriangle3Pts (const CVector& p0, const CVector& p1, const CVector& p2)
    {
        CVector cent;
        for (int i=0; i<N_DIM; ++i)
        {
            cent[i] = (p0[i] + p1[i] + p2[i]) / 3.;
        }
        
        return cent;
    }
    
    double Triangle::qualityScore (const vector<Point>& points, double aveTriArea, bool verbose, bool& passed)
    {
        // low score is better
        
        #include "Triangle.h"

        passed = true;
        
        //double a = mag( points[p[0]].dim - points[p[1]].dim );
        //double b = mag( points[p[0]].dim - points[p[2]].dim );
        //double c = mag( points[p[1]].dim - points[p[2]].dim );
        
        CVector p0 = points[p[0]].dim;
        CVector p1 = points[p[1]].dim;
        CVector p2 = points[p[2]].dim;
        
        double aveEdgeSize = sqrt ( (4./sqrt(3.) * aveTriArea) );
        CVector cnt;
        double radius;
        triPtsCircums (p0, p1, p2, cnt, radius);
        
        if (radius > (aveEdgeSize/2.))
        {
            passed = false;
        }
        
        /*double s = 0.5 * (a + b + c);        
        double r = (a*b*c) / (4. * sqrt( s*(s-a)*(s-b)*(s-c) ) );
        double areaEqui = 3*sqrt(3.)/4. * pow(r, 2.);
        double areaCurrent = areaTriangle (*this, points);
        
        // skewness
        double skew = 1. - areaCurrent/areaEqui;
        
        // area
        double areaSmall = min (areaCurrent, aveTriArea);
        double areaLarge = max (areaCurrent, aveTriArea);        
        double devAveArea = areaLarge / areaSmall;
        //double devAveArea = fabs(areaCurrent - aveTriArea) / aveTriArea;
        
        // aspect ratio
        double lon = max (max(a, b), c);
        double sho = min (min(a, b), c);
        //double aR = lon/sho;
        double aR = 1. - sho/lon; // min:0 (sho<<lon) || max:1 (sho=lon)
        //double aR = (lon-sho) / sho; // min:0 (lon=sho) || max:1 (lon=2*sho)
        
        double rSkew = cSkew*skew;
        //double rArea = cAA*devAveArea;
        double rAR = cAR*aR;
        //double score = rSkew + rArea + rAR;
        double score = rSkew + rAR;
        
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
            cout << "skew" << endl;
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
            cout << "ar" << endl;
            return score;
        }
        
        if (verbose)
        {
            cout << "skew = " << skew << endl;
            //cout << "devAveArea = " << devAveArea << endl;
            cout << "aR = " << aR << endl;
            cout << "rSkew = " << rSkew << endl;
            //cout << "rArea = " << rArea << endl;
            cout << "rAR = " << rAR << endl;
            cout << "score = " << score << endl;
        }
        
        return score;*/
    }
    
    bool triQuality (const CVector& p0, const CVector& p1, const CVector& p2, double rho)
    {
        CVector cnt;
        double radius;
        
        triPtsCircums (p0, p1, p2, cnt, radius);
        
        if (radius > rho)
        {
            return false;
        }
        
        return true;
    }
}


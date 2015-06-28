#include "AFT.h"

namespace AFT
{
    Edge::Edge()
    {
        newlyCreated = false;
    }
    
    void EdgeADT::build (const vector<Point>& points, const vector<Edge>& edges)
    {
        this->points.resize ( edges.size() );
        
        for (unsigned int e=0; e<edges.size(); ++e)
        {
            const Point& t0 = points[ edges[e].t[0] ];
            const Point& t1 = points[ edges[e].t[1] ];
            
            this->points[e] = this->createADTPoint (t0, t1);
        }
        
        ADT::build();
    }
    
    ADT::ADTPoint EdgeADT::createADTPoint (const Point& a, const Point& b)
    {
        ADTPoint vec;

        for (unsigned int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = min( a.dim[i], b.dim[i] );
            vec.dim[i*2+1] = max( a.dim[i], b.dim[i] );
        }

        vec.vertices.push_back (a);
        vec.vertices.push_back (b);

        return vec;
    }
    
    bool EdgeADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
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
    
        return true;
    }
    
    bool EdgeADT::compareFunction (const Node* node, const ADTPoint& targetPoint)
    {
        bool inside = doIntersect(targetPoint.vertices[0].dim, targetPoint.vertices[1].dim, node->p.vertices[0].dim, node->p.vertices[1].dim);

        return inside;
    }
    
    Edge createEdge (int indexA, int indexB, int belonging, bool newlyCreated)
    {
        Edge edge;
                
        edge.t.push_back (indexA);
        edge.t.push_back (indexB);
        edge.belonging = belonging;
        edge.newlyCreated = newlyCreated;

        return edge;
    }
    
    int edgeExists (const int ip1, const int ip2, const vector<Point>& points, const vector<Edge>& edges, bool& exists)
    {
        exists = false;
        
        for (unsigned int e=0; e<edges.size(); ++e)
        {
            int it0 = edges[e].t[0];
            int it1 = edges[e].t[1];
            
            if ( (it0 == ip1) || (it1 == ip1) )
            {
                if ( (it0 == ip2) || (it1 == ip2) )
                {
                    exists = true;
                    return e;
                }
            }
        }
    }
    
    bool checkEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT)
    {
        int result = -1;

        ADT::ADTPoint vec = edgeADT.createADTPoint (frontListPoint, closestPoint);

        result = edgeADT.search (vec);

        if (result != -1)
        {
            return true;
        }
        
        return false;
    }
    
    void knowParentTriangles (vector<Edge>& edges, const vector<Triangle>& triangles)
    {
        for (unsigned int t=0; t<triangles.size(); ++t)
        {
            for (int d=0; d<3; ++d)
            {
                edges[triangles[t].e[d]].nei.push_back (t);
            }
        }
        
        for (Edge& e: edges)
        {
            if (e.nei.size() == 1)
            {
                e.nei.push_back (-1);
            }
            else if (e.nei.size() == 0)
            {
                cout << "e.nei.size() == 0 in knowParentTriangles(...)" << endl;
                exit(-2);
            }
        }
    }
}

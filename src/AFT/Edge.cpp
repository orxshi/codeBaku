#include "AFT.h"

namespace AFT
{
    Edge::Edge()
    {
        newlyCreated = false;
        alive = true;
    }
    
    void EdgeADT::build (const vector<Point>& points, const vector<Edge>& edges)
    {
        this->points.resize ( edges.size() );
        
        for (unsigned int e=0; e<edges.size(); ++e)
        {
            const Point& t0 = points[ edges[e].t[0] ];
            const Point& t1 = points[ edges[e].t[1] ];
            
            this->points[e] = this->createADTPoint (t0, t1);
            this->points[e].idx = e;
        }
        
        ADT::build();
    }
    
    ADT::ADTPoint EdgeADT::createADTPoint (const AFT::Point& a, const AFT::Point& b)
    {
        ADTPoint vec;

        for (unsigned int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = min( a.dim[i], b.dim[i] );
            vec.dim[i*2+1] = max( a.dim[i], b.dim[i] );
        }

        vec.vertices.push_back (a.dim);
        vec.vertices.push_back (b.dim);

        return vec;
    }
    
    bool EdgeADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
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
    
        return true;
    }
    
    bool EdgeADT::compareFunction (const Node* node, const ADTPoint& targetPoint)
    {
        bool exactMatch;        
        bool inside = doIntersect (targetPoint.vertices[0], targetPoint.vertices[1], node->p->vertices[0], node->p->vertices[1], exactMatch);
        
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
    
    bool checkEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT, const vector<Edge>& edges, const vector<Point>& points, bool& exactMatch, int& result)
    {
        result = -1;

        ADT::ADTPoint vec = edgeADT.createADTPoint (frontListPoint, closestPoint);

        edgeADT.searchForNIntersections = false;
        result = edgeADT.search (vec);

        exactMatch = false;
        if (result != -1)
        {
            //cout << result << endl;
            //cout << edges.size() << endl;
            doIntersect (frontListPoint.dim, closestPoint.dim, points[edges[result].t[0]].dim, points[edges[result].t[1]].dim, exactMatch);
            return true;
        }
        
        return false;
    }
    
    int checkNumberOfEdgeIntersection (const Point& closestPoint, const Point& frontListPoint, EdgeADT& edgeADT)
    {
        ADT::ADTPoint vec = edgeADT.createADTPoint (frontListPoint, closestPoint);

        edgeADT.searchForNIntersections = true;
        edgeADT.search (vec);
        
        return edgeADT.nIntersections;
    }
    
    void knowParentTriangles (vector<Edge>& edges, const vector<Triangle>& triangles)
    {
        /*for (unsigned int t=0; t<triangles.size(); ++t)
        {
            for (int d=0; d<3; ++d)
            {
                edges[triangles[t].e[d]].nei.push_back (t);
            }
        }*/
        
        /*for (Edge& e: edges)
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
        }*/
    }
    
    void addToEdgeList (Edge& edge, int iP1, int iP2, vector<Edge>& edges, EdgeADT& edgeADT, vector<Point>& points)
    {
        bool tempBool;
        
        edges.push_back (edge);
        ADT::ADTPoint vec = edgeADT.createADTPoint (points[iP1], points[iP2]);
        vec.idx = edges.size() - 1;
        
        edgeADT.insert (vec, edgeADT.root, tempBool);
        
        points[iP1].e.push_back (edges.size() - 1);
        points[iP2].e.push_back (edges.size() - 1);
    }
    
    void eraseDeadEdges (vector<Edge>& edges, vector<Triangle>& triangles, vector<Point>& points)
    {
        for (int ie=0; ie<edges.size(); ++ie)
        {
            Edge& e = edges[ie];
            
            e.id = ie;
        }
        
        edges.erase (remove_if (edges.begin(), edges.end(), [](Edge& e) { return e.alive == false; }), edges.end());
        
        /*for (int ie=0; ie<edges.size(); ++ie)
        {
            Edge& e = edges[ie];
        
            if (e.alive == false)
            {
                edges.erase(edges.begin() + ie);
            }
        }*/
        
        for (int ie=0; ie<edges.size(); ++ie)
        {
            Edge& e = edges[ie];
            
            for (int ip: e.t)
            {
                Point& p = points[ip];
                
                for (int pp=0; pp<p.e.size(); ++pp)
                {
                    if (p.e[pp] == e.id)
                    {
                        p.e[pp] = ie;
                        break;
                    }
                }
            }
            
            for (int t: e.nei)
            {            
                if (t != -1)
                {
                    Triangle& tri = triangles[t];
                    
                    for (int pp=0; pp<tri.p.size(); ++pp)
                    {
                        if (tri.e[pp] == e.id)
                        {
                            tri.e[pp] = ie;
                            break;
                        }
                    }
                }
            }
        }
        
        for (int ie=0; ie<edges.size(); ++ie)
        {
            Edge& e = edges[ie];
            
            e.id = ie;
        }
    }
    
    bool Edge::eraseParentTri (int iTri)
    {
        for (int it=0; it<nei.size(); ++it)        
        {
            int t = nei[it];
        
            if (t == iTri)
            {
                nei.erase (nei.begin() + it);
                return true;
            }
        }
        
        return false;
    }
}

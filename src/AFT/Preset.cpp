#include "AFT.h"

namespace AFT
{
    void setPointsEdges (const vector<Grid>& gr, vector<Point>& points, vector<Edge>& edges, int newGridId)
    {
        vector<int> indices;
        int counter;

        for (int g=0; g<gr.size(); ++g)
        {
            for (const Face& f: gr[g].face)
            {
                //if (f.bouType == face_t::INTERIOR)
                //{
                    const Cell& L = gr[g].cell[f.nei[0]];
                    const Cell& R = gr[g].cell[f.nei[1]];

                    if (L.iBlank == iBlank_t::FIELD && R.iBlank == iBlank_t::FIELD)
                    {
                        if ( (L.trim == true && R.trim == false) || (R.trim == true && L.trim == false) )
                        {
                            // this face is quadrilateral because all grids have just one layer in z-direction

                            indices.clear();
                            for (const int iv: f.vtx)
                            {
                                Point v = gr[g].pt[iv];
                                //v.belonging = newGridId;

                                if (v.dim[2] == 0.)
                                {
                                    for (unsigned int p=0; p<points.size(); ++p)
                                    {
                                        counter = 0;

                                        for (int d=0; d<N_DIM; ++d)
                                        {
                                            if ( points[p].dim[d] == v.dim[d] )
                                            {
                                                ++counter;
                                            }
                                        }

                                        if (counter == N_DIM) // exists
                                        {
                                            indices.push_back (p);
                                            break;
                                        }
                                    }

                                    if (counter != N_DIM) // does not exist
                                    {
                                        points.push_back (v);
                                        indices.push_back ( points.size()-1 );
                                    }
                                }
                            }

                            Edge tempEdge = createEdge (indices[0], indices[1], g, false);
                            edges.push_back (tempEdge);
                        }
                    }
                //}
            }
        }
    }
    
    void createFrontList (const vector<Edge>& edges, vector<FrontMember>& frontList, const vector<Point>& points)
    {
        for (unsigned int e=0; e<edges.size(); ++e)
        {
            addToFrontList (e, frontList);
        }

        sortFrontList (frontList, points, edges);
        
        /*for (int e=0; e<edges.size(); ++e)
        {
            cout << edges[e].newlyCreated << endl;
        }
        
        exit(-2);*/
    }
}

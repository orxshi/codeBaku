#include "Cell.h"

Cell::Cell()
{
    vol    = 0.;
    iBlank = iBlank_t::UNDEFINED;
    trim   = false;
    newlyCreated = false;
    type = elmType_t::UNDEFINED;
    bc = BC::UNDEFINED;
    donor = NULL;
    receiver = NULL;
    for (Vector<N_DIM>& i: grad)
    {
        i.fill (0.);
    }
}

void Cell::set_centroid (const vector<Point>& pt)
{
    for (unsigned int i=0; i<cnt.size(); ++i)
    {
        cnt[i] = 0.;
    }

    for (const int p: vtx)
    {
        cnt += pt[p].dim;
    }

    cnt /= vtx.size();
}

void Cell::cons_to_prim()
{
    double k,ie;    
    
    prim[0] = cons[0];
    prim[1] = cons[1] / prim[0];
    prim[2] = cons[2] / prim[0];
    prim[3] = cons[3] / prim[0];

    k = 0.5 * ( pow(prim[1],2) + pow(prim[2],2) + pow(prim[3],2) );
    ie = cons[4] / prim[0] - k;

    prim[4] = prim[0] * (GAMMA - 1.) * ie;
}

void Cell::prim_to_cons()
{
    double k,ie;
    
    cons[0] = prim[0];
    cons[1] = prim[0] * prim[1];
    cons[2] = prim[0] * prim[2];
    cons[3] = prim[0] * prim[3];

    k = 0.5 * ( pow(prim[1],2) + pow(prim[2],2) + pow(prim[3],2) );
    ie = prim[4] / ( (GAMMA-1.) * prim[0] );

    cons[4] = prim[0] * (k + ie);
}

void Cell::interpolate()
{
    CVector dis;    
    
    if ( iBlank == iBlank_t::FRINGE )
    {
        for (int i=0; i<N_VAR; ++i)
        {
            dis = cnt - (*donor).cnt;

            prim[i] = (*donor).prim[i] + dotP((*donor).grad[i], dis);
        }

        prim_to_cons();
    }
}
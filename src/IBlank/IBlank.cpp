#include "../Grid/Grid.h"

void Grid::CellADT::build (const Grid& gr)
{
    //points.resize (gr.cell.size());
    points.resize (gr.n_in_elm);
    
    for (unsigned int c=0; c<points.size(); ++c)
    {
        for (int i=0; i<ADT_DIM; ++i)
        {
            points[c].dim[i*2]   = BIG_POS_NUM;
            points[c].dim[i*2+1] = BIG_NEG_NUM;
        }

        for (const int ip: gr.cell[c+gr.n_bou_elm].vtx)
        {
            const Point& p = gr.pt[ip];
            
            for (int i=0; i<ADT_DIM; ++i)
            {
                points[c].dim[i*2]   = min (p.dim[i], points[c].dim[i*2]);
                points[c].dim[i*2+1] = max (p.dim[i], points[c].dim[i*2+1]);
            }
            
            /*points[c].dim[0] = min (p.dim[0], points[c].dim[0]);
            points[c].dim[1] = max (p.dim[0], points[c].dim[1]);

            points[c].dim[2] = min (p.dim[1], points[c].dim[2]);
            points[c].dim[3] = max (p.dim[1], points[c].dim[3]);

            points[c].dim[4] = min (p.dim[2], points[c].dim[4]);
            points[c].dim[5] = max (p.dim[2], points[c].dim[5]);*/

            points[c].vertices.push_back (p);
        }

        points[c].idx = c+gr.n_bou_elm;
    }
    
    ADT::build();
}

void Grid::setWallDistance (const Solid& sld)
{
    for (int c=0; c<cell.size(); ++c)
    {
        double d = BIG_POS_NUM;
        
        for (int g=0; g<n_bou_elm; ++g)
        {
            if (cell[g].phys == sld.phys)
            {
                CVector tmpv = cell[g].cnt - cell[c].cnt;

                d = min ( d, mag (tmpv) );
            }
        }

        cell[c].wallDistance = d;
    }
}

void Grid::identifyIBlank (Grid& gr)
{
    function<int(Cell&)> getIndex = [&] (Cell& cll)
    {
        ADT::ADTPoint vec;
        CVector cnt;
        
        cnt = cll.cnt;
        
        for (int i=0; i<ADT_DIM; ++i)
        {
            vec.dim[i*2]   = cnt[i];
            vec.dim[i*2+1] = cnt[i];
        }

        /*vec.dim[0] = cnt[0];
        vec.dim[2] = cnt[1];
        vec.dim[4] = cnt[2];

        vec.dim[1] = cnt[0];
        vec.dim[3] = cnt[1];
        vec.dim[5] = cnt[2];*/
        
        return gr.cellADT.search (vec);
    };
            
    for (int c=n_bou_elm; c<cell.size(); ++c)
    {
        int index;
        
        if (cell[c].iBlank == iBlank_t::UNDEFINED)
        {            
            index = getIndex (cell[c]);
            
            if (index == -1)
            {                
                cell[c].iBlank = iBlank_t::FIELD;
            }
            else
            {             
                cell[c].iBlank = iBlank_t::FRINGE;
                
                if (gr.cell[index].wallDistance < cell[c].wallDistance)
                {
                    cell[c].iBlank = iBlank_t::FRINGE;
                    gr.cell[index].iBlank = iBlank_t::FIELD;
                    cell[c].donor = &gr.cell[index];
                    gr.cell[index].receiver = &cell[c];
                }
                else
                {
                    cell[c].iBlank = iBlank_t::FIELD;
                    gr.cell[index].iBlank = iBlank_t::FRINGE;
                    gr.cell[index].donor = &cell[c];
                    cell[c].receiver = &gr.cell[index];
                }
            }
        }
    }
    
    for (int c=0; c<n_bou_elm; ++c)
    {
        int index;
        
        if (cell[c].iBlank == iBlank_t::UNDEFINED)
        {
            index = getIndex (cell[c]);
        
            if (index == -1)
            {
                cell[c].iBlank == iBlank_t::NA;
            }
            else
            {
                cell[c].iBlank = iBlank_t::FRINGE;
                cell[c].donor = &gr.cell[index];
                gr.cell[index].receiver = &cell[c];
                gr.cell[index].iBlank = iBlank_t::FIELD;
            }
        }
    }
}

bool Grid::CellADT::compareFunction (const Node* node, const ADTPoint& targetPoint)
{
    bool inside;
    CVector tempVector;
    unsigned int size = node->p.vertices.size();
    double frac[size];
    vector<CVector> xc(size);

    for (unsigned int v=0; v<size; ++v)
    {
        xc[v] = node->p.vertices[v].dim;
    }
    
    for (int i=0; i<ADT_DIM; ++i)
    {
        tempVector[i] = targetPoint.dim[i*2];
    }
    
    /*tempVector[0] = targetPoint.dim[0];
    tempVector[1] = targetPoint.dim[2];
    tempVector[2] = targetPoint.dim[4];*/
    
    osInterpolants (xc, tempVector, size, frac);

    inside = true;
    for (unsigned int v=0; v<size; ++v)
    {
        if (frac[v]<=0 || frac[v]>=1)
        {
            inside = false;
            break;
        }
    }

    return inside;
}

bool Grid::CellADT::doCubesOverlap (const Node* node, const ADTPoint& targetPoint)
{
    bool insideCube = true;

    for (int d=0; d<ADT_DIM; ++d)
    {
        if (!(node->p.dim[d*2] <= targetPoint.dim[d*2]) || !(node->p.dim[d*2+1] >= targetPoint.dim[d*2+1]))
        {
            insideCube = false;
            break;
        }
    }

    return insideCube;
}

void Grid::interpolate()
{
    for (Cell& cll: cell) { cll.interpolate(); }
}
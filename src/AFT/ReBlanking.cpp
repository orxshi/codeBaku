#include "AFT.h"

void Grid::fieldToFringe (int crt)
{
    int counter;

    for (Cell& c: cell)
    {
        if (c.iBlank == iBlank_t::FIELD)
        {
            counter = 0;
            for (const int nei: c.nei)
            {
                Cell& n = cell[nei];
                
                if (n.iBlank == iBlank_t::FRINGE)
                {
                    ++counter;
                    if (counter == crt)
                    {
                        c.receiver->iBlank = iBlank_t::FIELD;
                        c.receiver->donor = NULL;
                        c.receiver->receiver = &c;

                        c.iBlank = iBlank_t::FRINGE;
                        c.donor = c.receiver;
                        c.receiver = NULL;
                        break;
                    }
                }
            }
        }
    }
}
    
void Grid::fringeToField (int crt)
{
    int counter;
    int donor;

    for (Cell& c: cell)
    {
        if (c.iBlank == iBlank_t::FRINGE)
        {
            counter = 0;            
            for (const int nei: c.nei)
            {
                Cell& n = cell[nei];
                
                if (n.iBlank == iBlank_t::FIELD)
                {
                    ++counter;
                    if (counter == crt)
                    {
                        c.donor->iBlank = iBlank_t::FRINGE;
                        c.donor->donor = &c;
                        c.donor->receiver = NULL;

                        c.iBlank = iBlank_t::FIELD;
                        c.receiver = c.donor;
                        c.donor = NULL;

                        break;
                    }
                }
            }
        }
    }
}

void Grid::blankWithPhys (int phys)
{
    for (int c=0; c<n_bou_elm; ++c)
    {
        Cell& cll = cell[c];
        
        if (cll.phys == phys)
        {
            cll.iBlank = iBlank_t::FRINGE;
            cell[cll.nei[0]].iBlank = iBlank_t::FRINGE;
        }
    }
}
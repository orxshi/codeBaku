#include "ADT.h"

ADT::ADT()
{
    sizeTree = 0;
    root = NULL;
    nIntersections = 0;
    searchForNIntersections = false;
}

ADT::~ADT()
{
    destroy_tree();
}

void ADT::destroy_tree(Node *&leaf)
{
    if (leaf != NULL)
    {
        destroy_tree(leaf->left);
        destroy_tree(leaf->right);
        delete leaf;
        leaf = NULL;
    }
}

void ADT::destroy_tree()
{
    destroy_tree(root);

    for (unsigned int i=0; i<searchStack.size(); ++i)
    {
        searchStack[i] = NULL;
    }
}

void ADT::insert (const ADTPoint& point, Node* node, bool& isInserted)
{
    unsigned int j = node->level % ADT_VAR;
    node->key = 0.5 * (node->c[j] + node->d[j]);

    // point is in left region
    if (point.dim[j] < node->key)
    {
        if (node->left == NULL)
        {
            node->left = new Node();
            ++sizeTree;

            for (int d=0; d<ADT_VAR; ++d)
            {
                if (d != j)
                {
                    node->left->c[d] = node->c[d];
                    node->left->d[d] = node->d[d];
                }
                else
                {
                    node->left->c[d] = node->c[d];
                    node->left->d[d] = node->key;
                }
            }

            node->left->level = node->level + 1;
            isInserted = true;
            node->left->p = point;
        }
        else
        {
            insert(point, node->left, isInserted);
        }
    }
    // point is in right region
    else
    {
        if (node->right == NULL)
        {
            node->right = new Node();
            ++sizeTree;

            for (int d=0; d<ADT_VAR; ++d)
            {
                if (d != j)
                {
                    node->right->c[d] = node->c[d];
                    node->right->d[d] = node->d[d];
                }
                else
                {
                    node->right->c[d] = node->key;
                    node->right->d[d] = node->d[d];
                }
            }

            node->right->level = node->level + 1;
            isInserted = true;
            node->right->p = point;
        }
        else
        {
            insert(point, node->right, isInserted);
        }
    }
}

bool ADT::doRegionOverlap (const unsigned int j, const Node* node, const ADTPoint& targetPoint)
{
    bool regionOverlap = true;

    if ( (j % 2 == 0) )
    {
        if ( !(node->c[j] <= targetPoint.dim[j+1]) )
        {
            regionOverlap = false;
        }
    }
    else
    {
        if ( !(node->d[j] >= targetPoint.dim[j-1]) )
        {
            regionOverlap = false;
        }
    }

    return regionOverlap;
}

void ADT::searchChildren (const Node* node, const ADTPoint& targetPoint)
{
    unsigned int j = node->level % ADT_VAR;

    // check left
    if (node->left != NULL)
    {
        if ( doRegionOverlap (j, node->left, targetPoint) )
        {
            searchStack.push_back (node->left);
        }
    }

    // check right
    if (node->right != NULL)
    {
        if ( doRegionOverlap (j, node->right, targetPoint) )
        {
            searchStack.push_back (node->right);
        }
    }
}

void ADT::search (Node* node, const ADTPoint& targetPoint, int& index)
{
    if (node != NULL)
    {
        //cout << "node->p.dim[0] = " << node->p.dim[0] << endl;
        //cout << "node->p.dim[1] = " << node->p.dim[1] << endl;
        
        //cout << "targetPoint.dim[0] = " << targetPoint.dim[0] << endl;
        //cout << "targetPoint.dim[1] = " << targetPoint.dim[1] << endl;
        
        // check whether the point is inside the element
        if ((node->p.idx!=-1) && doCubesOverlap (node, targetPoint) && compareFunction (node, targetPoint) )
        {
            if (searchForNIntersections)
            {
                ++nIntersections;
                ids.push_back (node->p.idx);
                
                searchChildren (node, targetPoint);

                if (!searchStack.empty())
                {
                    Node* last = searchStack.back();                
                    searchStack.back() = NULL;                
                    searchStack.pop_back();
                    search (last, targetPoint, index);
                }
            }
            else
            {
                fill (searchStack.begin(), searchStack.end(), nullptr);
                searchStack.clear();
            
                index = node->p.idx;
                //cout << "node->p.idx = " << node->p.idx << endl;
                //cout << "node->p.dim[0] = " << node->p.dim[0] << endl;
                //cout << "node->p.dim[1] = " << node->p.dim[1] << endl;
            
                if (node != root)
                {
                    node = NULL;
                }
            }
        }
        else
        {
            searchChildren (node, targetPoint);

            if (!searchStack.empty())
            {
                Node* last = searchStack.back();                
                searchStack.back() = NULL;                
                searchStack.pop_back();
                search (last, targetPoint, index);
            }
        }
    }
}

int ADT::search (const ADTPoint& targetPoint)
{
    int i = -1;
    nIntersections = 0;
    ids.clear();
    
    fill (searchStack.begin(), searchStack.end(), nullptr);
    searchStack.clear();

    bool regionOverlap = true;

    for (int d=0; d<ADT_DIM; ++d)
    {
        if ( !(root->c[2*d] <= targetPoint.dim[2*d+1]) || !(root->d[2*d+1] >= targetPoint.dim[2*d]) )
        {
            regionOverlap = false;
        }
    }

    if (regionOverlap)
    {
        search (root, targetPoint, i);
    }
    /*else
    {
        cout << "no overlap" << endl;
        cout << "root->c[0] = " << root->c[0] << endl;
        cout << "root->c[1] = " << root->c[1] << endl;
        cout << "root->c[2] = " << root->c[2] << endl;
        cout << "root->c[3] = " << root->c[3] << endl;
        cout << "root->c[4] = " << root->c[4] << endl;
        cout << "root->c[5] = " << root->c[5] << endl;
        
        cout << "root->d[0] = " << root->d[0] << endl;
        cout << "root->d[1] = " << root->d[1] << endl;
        cout << "root->d[2] = " << root->d[2] << endl;
        cout << "root->d[3] = " << root->d[3] << endl;
        cout << "root->d[4] = " << root->d[4] << endl;
        cout << "root->d[5] = " << root->d[5] << endl;
        
        cout << "targetPoint.dim[0] = " << targetPoint.dim[0] << endl;
        cout << "targetPoint.dim[1] = " << targetPoint.dim[1] << endl;
        cout << "targetPoint.dim[2] = " << targetPoint.dim[2] << endl;
        cout << "targetPoint.dim[2] = " << targetPoint.dim[3] << endl;
        cout << "targetPoint.dim[2] = " << targetPoint.dim[4] << endl;
        cout << "targetPoint.dim[2] = " << targetPoint.dim[5] << endl;
    }*/

    fill (searchStack.begin(), searchStack.end(), nullptr);
    searchStack.clear();    
    
    return i;
}

void ADT::build ()
{
    bool isInserted;
    
    if (root == NULL)
    {
        root = new Node();
        ++sizeTree;
        root->level = 0;

        for (int d=0; d<ADT_VAR; ++d)
        {
            root->c[d] = BIG_POS_NUM;
            root->d[d] = BIG_NEG_NUM;
        }

        for (unsigned int p=0; p<points.size(); ++p)
        {
            for (int d=0; d<ADT_VAR; ++d)
            {
                root->c[d] = min (points[p].dim[d], root->c[d]);
                root->d[d] = max (points[p].dim[d], root->d[d]);
            }
        }

        root->p = points.front();

        points.erase (points.begin());
    }

    while (!points.empty())
    {
        isInserted = false;
        insert (points.front(), root, isInserted);
        if (isInserted)
        {
            points.erase (points.begin());
        }
    }
}
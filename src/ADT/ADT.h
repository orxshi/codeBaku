/* 
 * File:   ADT.h
 * Author: Orhan Shibliyev
 *
 * Created on July 16, 2014, 6:32 PM
 */

#ifndef ADT_H
#define	ADT_H

#include <vector>
#include <algorithm>
#include <iostream>
#include "../Vector/Vector.h"
#include "../Constants.h"
#include "../Point/Point.h"

#define ADT_DIM 3
#define ADT_VAR (2*ADT_DIM)

using std::vector;
using std::min;
using std::max;
using std::fill;
using std::cout;
using std::endl;
using std::cin;

struct ADT
{
    struct ADTPoint
    {
        int idx;
        vector <Point> vertices;
        Vector <ADT_VAR> dim;
        /*
         * 0: xmin
         * 1: xmax
         * 2: ymin
         * 3: ymax
         * 4: zmin
         * 5: zmax
         */
    };
    
    unsigned int sizeTree;
    vector<ADTPoint> points;    

    struct Node
    {
        Vector<ADT_VAR> c, d;
        ADTPoint p;
        double key;
        Node* left;
        Node* right;
        unsigned int level;

        Node()
        {
            left = NULL;
            right = NULL;
        }
    };

    Node *root;
    vector <Node*> searchStack;    

    void destroy_tree (Node *&leaf);
    void destroy_tree();

    void insert (const ADTPoint& point, Node *node, bool& isInserted);

    bool doRegionOverlap (const unsigned int j, const Node *node, const ADTPoint& targetPoint);

    virtual bool compareFunction (const Node *node, const ADTPoint& targetPoint) {}
    virtual bool doCubesOverlap (const Node *node, const ADTPoint& targetPoint) {}

    void searchChildren (const Node *node, const ADTPoint& targetPoint);

    void search (Node* node, const ADTPoint& targetPoint, int& index);
    int search (const ADTPoint& targetPoint);

    // Constructor
    ADT();
    
    // Destructor
    ~ADT();

    void build();
};

#endif	/* ADT_H */


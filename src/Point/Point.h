/* 
 * File:   Point.h
 * Author: orxan
 *
 * Created on August 26, 2014, 8:36 PM
 */

#ifndef POINT_H
#define	POINT_H

#include "../Vector/Vector.h"

struct Point
{
    CVector dim;
    int belonging;
    bool newlyCreated;
    
    Point();
};

#endif	/* POINT_H */


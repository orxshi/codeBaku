#include "AFT.h"

namespace AFT
{
    bool doIntersect (const CVector& p1, const CVector& q1, const CVector& p2, const CVector& q2)
    {
        // http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/

        // The main function that returns true if line segment 'p1q1'
        // and 'p2q2' intersect

        // Find the four orientations needed for general and
        // special cases

        int o1 = orientation(p1, q1, p2);
        int o2 = orientation(p1, q1, q2);
        int o3 = orientation(p2, q2, p1);
        int o4 = orientation(p2, q2, q1);    

        if (o1 == -1 || o2 == -1 || o3 == -1 || o4 == -1)
        {
            return true;
        }
        // General case
        else if (o1 != o2 && o3 != o4)
        {
            if (o1*o2*o3*o4 != 0)
            {
                return true;
            }
        }
        else if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0)
        {
            // Special Cases
            // p1, q1 and p2 are colinear and p2 lies on segment p1q1
            if (o1 == 0 && onSegment(p1, p2, q1)) return true;

            // p1, q1 and p2 are colinear and q2 lies on segment p1q1
            if (o2 == 0 && onSegment(p1, q2, q1)) return true;

            // p2, q2 and p1 are colinear and p1 lies on segment p2q2
            if (o3 == 0 && onSegment(p2, p1, q2)) return true;

             // p2, q2 and q1 are colinear and q1 lies on segment p2q2
            if (o4 == 0 && onSegment(p2, q1, q2)) return true;

            unsigned int count = 0;

            if ( onSegment(p1, p2, q1) ) ++count;
            if ( onSegment(p1, q2, q1) ) ++count;
            if ( onSegment(p2, p1, q2) ) ++count;
            if ( onSegment(p2, q1, q2) ) ++count;

            if (count == 0)
            {
                if ( ((p1[0] == p2[0]) && (p1[1] == p2[1])) || ((p1[0] == q2[0]) && (p1[1] == q2[1])) )
                {
                    if ( ((q1[0] == p2[0]) && (q1[1] == p2[1])) || ((q1[0] == q2[0]) && (q1[1] == q2[1])) )
                    {
                        return true;
                    }
                }
            }
        }

        return false; // Doesn't fall in any of the above cases
    }
    
    int orientation (const CVector& p, const CVector& q, const CVector& r)
    {
        // To find orientation of ordered triplet (p, q, r).
        // The function returns following values
        // 0 --> p, q and r are colinear
        // 1 --> Clockwise
        // 2 --> Counterclockwise

        // See 10th slides from following link for derivation of the formula
        // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf

        double val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);

        if (val == 0.) return 0;  // colinear

        return (val > 0.) ? 1 : 2; // clock or counterclock wise
    }
    
    bool onSegment (const CVector& p, const CVector& q, const CVector& r)
    {
        // Given three colinear Vectors p, q, r, the function checks if
        // Vector q lies on line segment 'pr'

        bool AllEqualX;
        bool AllEqualY;

        if ( q[0] == max(p[0], r[0]) && q[0] == min(p[0], r[0]) )
        {
            AllEqualX = true;
        }
        else
        {
            AllEqualX = false;
        }

        if ( q[1] == max(p[1], r[1]) && q[1] == min(p[1], r[1]) )
        {
            AllEqualY = true;
        }
        else
        {
            AllEqualY = false;
        }

        if ( (((q[0] < max(p[0], r[0])) && (q[0] > min(p[0], r[0]))) || AllEqualX) &&
             (((q[1] < max(p[1], r[1])) && (q[1] > min(p[1], r[1]))) || AllEqualY) )
           return true;

        return false;
    }
}

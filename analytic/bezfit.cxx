#include <stdio.h>
#include <dlib/optimization.h"

#define _DMAX 1000000

typedef dlib::matrix<double,0,1> column_vector;

struct Point
{
    double x;
    double y;
};

static Point points[] =
{
    { 254 , 102 },
    { 226 , 63  },
    { 185 , 49  },
    { 146 , 74  },
    { 142 , 119 },
    { 117 , 169 },
    { 86  , 214 },
    { 40  , 200 },
};

//
// cubicBezier
//
// p1 - start point
// c1 - first control point
// c2 - second control point
// p2 - end point
//
double cubicBezier(double p1, double c1, double c2, double p2, double t)
{
    double s = (1 - t);

    double v = 0;
    v = v + (1 * p1 * s * s * s);
    v = v + (3 * c1 * s * s * t);
    v = v + (3 * c2 * s * t * t);
    v = v + (1 * p2 * t * t * t);

    return v;        
};

Point cubicBezier(double p1x, double p1y, double c1x, double c1y,
                  double c2x, double c2y, double p2x, double p2y,
                  double t)
{
    Point pt;
    pt.x = cubicBezier(p1x, c1x, c2x, p2x, t);
    pt.y = cubicBezier(p1y, c1y, c2y, p2y, t);
    return pt;
}

// Any distance function can be used for optimisation.  This one, where we want
// to find the least-squares is most common. 
double dist(double x1, double y1, double x2, double y2)
{
    double x = x2 - x1;
    double y = y2 - y1;

    return (x * x) + (y * y);
}

// This is the function that the optimiser calls repeatedly with different
// parameters, attempting to get the lowest return value it can.
double rateCurve(const column_vector& params)
{
    double p1x = points[0].x;
    double p1y = points[0].y;
    double c1x = params(0,0);
    double c1y = params(1,0);
    double c2x = params(2,0);
    double c2y = params(3,0);
    double p2x = points[7].x;
    double p2y = points[7].y;

    double distances = 0;

    for (Point target : points)
    {
        double distance = _DMAX;

        for (double t = 0; t <= 1; t += 0.02)
        {
            Point pt = cubicBezier( p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y, t);

            double dCandidate = dist(pt.x, pt.y, target.x, target.y);

            distance = std::min(distance, dCandidate);
        }

        distances += distance;
    }

    // Thats the curve-fitting done.  Now incentivise slightly smoother curves.
    double p1c1 = dist(p1x, p1y, c1x, c1y);
    double p2c2 = dist(p2x, p2y, c2x, c2y);

    return distances + pow(p1c1, 0.6) + pow(p2c2, 0.6);
}

int main(int argc, char* argv[])
{
    column_vector params(4);
    params = points[7].x, points[7].y, points[0].x, points[0].y;

    dlib::find_min_using_approximate_derivatives(
            dlib::cg_search_strategy(),
            dlib::objective_delta_stop_strategy(1).be_verbose(),
            rateCurve,
            params,
            -1);

    printf("p1x = %f;\n", points[0].x);
    printf("p1y = %f;\n", points[0].y);
    printf("c1x = %f;\n", params(0,0));
    printf("c1y = %f;\n", params(1,0));
    printf("c2x = %f;\n", params(2,0));
    printf("c2y = %f;\n", params(3,0));
    printf("p2x = %f;\n", points[7].x);
    printf("p2y = %f;\n", points[7].y);

    return 0;
}
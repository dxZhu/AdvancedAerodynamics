/*
 * FindRoot.cpp
 *
 *  Created on: May 12, 2017
 *      Author: dxzhu
 */

//using this part to get the last point on the C- line from wall

#include "main.h"

using namespace std;

void FindRoot(double xp, double yp, double& x3, double& y3, double& u3, double& v3)
{
    double dx = 1.3;
    double Md = 3.13;
    double alphap = asin(1/Md);
    x3 = xp + dx;
    y3 = yp + tan(alphap)*(dx);
    u3 = 1.9945;
    v3 = 0.00001;
}



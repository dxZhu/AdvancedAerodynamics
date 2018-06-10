/*
 * main.cpp
 *
 *  Created on: May 2, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

int main()
{
    int timestep = 500;
    int spacestep = 51;
    double dx;
    double dt;
    double gamma;
    double x0In, x0Out;
    double rho0In, rho0Out;
    double p0In, p0Out;
    double v0In, v0Out;
    double T0In, T0Out;

    arma::rowvec space;
    arma::mat velocity;
    arma::mat pressure;
    arma::mat density;
    arma::mat temperature;
    arma::rowvec area;

    //Initial Value for inlet and outlet
    x0In = 0.0;
    x0Out = 1.0;
    v0In = 0.183697;
    v0Out = 2.0;
    rho0In = 0.985999;
    rho0Out = 0.01;
    p0In = 0.980454;
    p0Out = 0.01;
    T0In = p0In/rho0In;
    T0Out = p0Out/rho0Out;
    dx = 1.0/(spacestep-1);
    dt = 0.005;
    //Flow properties
    gamma = 1.4;

    GetInitialValuePlane(x0In, x0Out, rho0In, rho0Out, p0In, p0Out, v0In, v0Out, T0In, T0Out, dx, spacestep, timestep,
		area, space, velocity, pressure, density, temperature);
for(int n=0; n<499; n++)
{
    GetPredictor(dx, dt, spacestep, n, gamma,
		area, velocity, pressure, density, temperature);

    GetCorrector(dx, dt, spacestep, n, gamma,
		area, velocity, pressure, density, temperature);
};

    PrintSolutions(area, space, velocity, pressure, density, temperature);

    ExtractSolutions(area, space, velocity, pressure, density, temperature);
}



/*
 * GetInitialValuePlane.cpp
 *
 *  Created on: May 2, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

void GetInitialValuePlane(double x0In, double x0Out, double rho0In, double rho0Out, double p0In, double p0Out, double v0In, double v0Out, double T0In, double T0Out,
	double dx, int spacestep, int timestep,
	arma::rowvec& area, arma::rowvec& space, arma::mat& velocity, arma::mat& pressure, arma::mat& density, arma::mat& temperature)
{

    //Initialize matrix
    space = arma::zeros<arma::rowvec>(spacestep);
    velocity = arma::zeros<arma::mat>(timestep, spacestep);
    pressure = arma::zeros<arma::mat>(timestep, spacestep);
    density = arma::zeros<arma::mat>(timestep, spacestep);
    temperature = arma::zeros<arma::mat>(timestep, spacestep);
    area = arma::zeros<arma::rowvec>(spacestep);

    int row = 0; //first time step
    double dx0 = (x0Out - x0In)/(spacestep-1);
    double dv0 = (v0Out - v0In)/(spacestep-1);
    double dp0 = (p0Out - p0In)/(spacestep-1);
    double drho0 = (rho0Out - rho0In)/(spacestep-1);
    double dT0 = (T0Out - T0In)/(spacestep-1);

    for(int i=0; i<spacestep; i++){
	space(i) = 0 + i*dx0;
	velocity(row, i) = v0In + i*dv0;
	pressure(row, i) = p0In + i*dp0;
	density(row, i) = rho0In + i*drho0;
	temperature(row, i) = T0In + i*dT0;
	area(i) = 1 + 10*(space(i)-0.5)*(space(i)-0.5);
    }


}



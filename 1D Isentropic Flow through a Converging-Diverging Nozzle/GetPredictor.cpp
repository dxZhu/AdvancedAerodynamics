/*
 * GetPredictor.cpp
 *
 *  Created on: May 2, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

void GetPredictor(double dx, double dt, int spacestep, int n, double gamma,
	arma::rowvec area, arma::mat& velocity, arma::mat& pressure, arma::mat& density, arma::mat& temperature)
{

    double dpdt;
    double dvdt;
    double drhodt;

    for(int i=0; i<spacestep-1; i++)
    {
	//change dt for each time step
	dpdt = -(velocity(n,i)*(pressure(n,i+1)-pressure(n,i))/dx + gamma*pressure(n,i)*(velocity(n,i+1)-velocity(n,i))/dx  +  gamma*velocity(n,i)*pressure(n,i)/area(i)*(area(i+1)-area(i))/dx);
	drhodt = -(velocity(n,i)*(density(n,i+1)-density(n,i))/dx + density(n,i)*(velocity(n,i+1)-velocity(n,i))/dx + density(n,i)*velocity(n,i)/area(i)*(area(i+1)-area(i))/dx);
	dvdt = -((gamma+1)/(2*gamma)/density(n,i)*(pressure(n,i+1)-pressure(n,i))/dx  +  velocity(n,i)*(velocity(n,i+1)-velocity(n,i))/dx);
	pressure(n+1,i) = pressure(n,i) + dt * dpdt;
	density(n+1,i) = density(n,i) + dt * drhodt;
	velocity(n+1,i) = velocity(n,i) + dt * dvdt;
	temperature(n+1,i) = pressure(n+1,i)/density(n+1,i);
    }

    pressure(n+1,0) = pressure(n,0);
    density(n+1,0) = density(n,0);
    velocity(n+1,0) = velocity(n,0);
    temperature(n+1,0) = temperature(n,0);
    pressure(n+1,spacestep-1) = pressure(n+1,spacestep-2) + (pressure(n+1,spacestep-2)-pressure(n+1,spacestep-3));
    density(n+1,spacestep-1) = density(n+1,spacestep-2) + (density(n+1,spacestep-2)-density(n+1,spacestep-3));
    velocity(n+1,spacestep-1) = velocity(n+1,spacestep-2) + (velocity(n+1,spacestep-2)-velocity(n+1,spacestep-3));
    temperature(n+1,spacestep-1) = temperature(n+1,spacestep-2) + (temperature(n+1,spacestep-2)-temperature(n+1,spacestep-3));

}



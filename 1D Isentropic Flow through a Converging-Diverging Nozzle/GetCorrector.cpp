/*
 * GetCorrector.cpp
 *
 *  Created on: May 2, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

void GetCorrector(double dx, double dt, int spacestep, int n, double gamma,
	arma::rowvec area, arma::mat& velocity, arma::mat& pressure, arma::mat& density, arma::mat& temperature)
{
    double dpdt;
    double drhodt;
    double dvdt;
    int Ni = 1000;
    double ep = 0.0001;
    arma::rowvec Tpressure;
    arma::rowvec Tdensity;
    arma::rowvec Ttemperature;
    arma::rowvec Tvelocity;
cout<< "timestep is " << n<< endl;
for(int j=0; j<1; j++){
    for(int i=1; i< spacestep; i++)
    {
	Tpressure = pressure.row(n+1);
	Tdensity = density.row(n+1);
	Tvelocity = velocity.row(n+1);
	Ttemperature = temperature.row(n+1);

	dpdt = -(velocity(n+1,i)*(pressure(n+1,i)-pressure(n+1,i-1))/dx  +  gamma*pressure(n+1,i)*(velocity(n+1,i)-velocity(n+1,i-1))/dx + gamma*velocity(n+1,i)*pressure(n+1,i)/area(i)*(area(i)-area(i-1))/dx);
	drhodt = -(velocity(n+1,i)*(density(n+1,i)-density(n+1,i-1))/dx  +  density(n+1,i)*(velocity(n+1,i)-velocity(n+1,i-1))/dx  +  density(n+1,i)*velocity(n+1,i)/area(i)*(area(i)-area(i-1))/dx);
	dvdt = -((gamma+1)/(2.0*gamma)/density(n+1,i)*(pressure(n+1,i)-pressure(n+1,i-1))/dx  +  velocity(n+1,i)*(velocity(n+1,i)-velocity(n+1,i-1))/dx);
	pressure(n+1, i) = 0.5*(pressure(n,i) + pressure(n+1,i) + dt*dpdt);
	density(n+1, i) = 0.5*(density(n,i) + density(n+1,i) + dt*drhodt);
	velocity(n+1, i) = 0.5*(velocity(n,i) + velocity(n+1,i) + dt*dvdt);
	temperature(n+1, i) = pressure(n+1,i)/density(n+1,i);
	}
    pressure(n+1,spacestep-1) = pressure(n+1,spacestep-2) + (pressure(n+1,spacestep-2)-pressure(n+1,spacestep-3));
    density(n+1,spacestep-1) = density(n+1,spacestep-2) + (density(n+1,spacestep-2)-density(n+1,spacestep-3));
    velocity(n+1,spacestep-1) = velocity(n+1,spacestep-2) + (velocity(n+1,spacestep-2)-velocity(n+1,spacestep-3));
    temperature(n+1,spacestep-1) = temperature(n+1,spacestep-2) + (temperature(n+1,spacestep-2)-temperature(n+1,spacestep-3));

    //Check for converge:
   /* if(abs((Tpressure(spacestep-1)-pressure(n+1,spacestep-1))/pressure(n+1,spacestep-1)) <ep)
    {break;}*/

}


}



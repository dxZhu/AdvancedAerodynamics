/*
 * Properties.cpp
 *
 *  Created on: Mar 10, 2017
 *      Author: dxzhu
 */
#include "main.h"

using namespace std;

//Get Pressure, StaticPressure, Density, temperature and Cp at cone surface for N discretized Wave angles:
void Properties(double Mach1, double gamma, double Vr, double Vpsi, double StaticPressure,
	double& Mach, double& Pressure, double& Temperature, double& Density, double& Cp)
{
    double MachStar;

    MachStar = sqrt(Vr*Vr+Vpsi*Vpsi);
    Mach = sqrt(2*MachStar*MachStar/(gamma+1)/(1 - (gamma-1)/(gamma+1)*MachStar*MachStar));
    Pressure = pow((1+(gamma-1)/2*Mach1*Mach1)/(1+(gamma-1)/2*Mach*Mach), gamma/(gamma-1)) * StaticPressure;
    Temperature = (1+(gamma-1)/2*Mach1*Mach1)/(1+(gamma-1)/2*Mach*Mach);
    Density = Pressure/Temperature;
    Cp = 2*(Pressure-1)/gamma/Mach1/Mach1;
}


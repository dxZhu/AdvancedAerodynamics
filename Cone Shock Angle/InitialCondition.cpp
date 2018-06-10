/*
 * InitialCondition.cpp
 *
 *  Created on: Mar 10, 2017
 *      Author: dxzhu
 */
#include "main.h"

using namespace std;

void InitialCondition(double WaveAngle, double Mach1, double Mach1Star, double Nu, double Ni, double gamma,
	double& dPsi, double& VelocityR1, double& VelocityPsi1, double& Temperature, double& Density, double& Pressure, double& StaticPressure)
{
    double Mach2Star; //M2*
    double Mach2;
    dPsi = -WaveAngle/Ni;
    VelocityR1 = Mach1Star * cos(WaveAngle);
    VelocityPsi1 = -Mach1Star * sin(WaveAngle) * (1 - 2.0/(gamma+1) * (Mach1*Mach1*sin(WaveAngle)*sin(WaveAngle)-1)/(Mach1*Mach1*sin(WaveAngle)*sin(WaveAngle)));
    Pressure = 1 + (2.0*gamma/(gamma+1))*(Mach1*Mach1*sin(WaveAngle)*sin(WaveAngle) - 1);
    Density = 1/(1 - 2.0/(gamma+1)*(Mach1*Mach1*sin(WaveAngle)*sin(WaveAngle) - 1)/(Mach1*Mach1*sin(WaveAngle)*sin(WaveAngle)));
    Temperature = 1 + 2.0*(gamma-1)/(gamma+1)/(gamma+1)*(Mach1*Mach1*sin(WaveAngle)*sin(WaveAngle)-1)/(Mach1*Mach1*sin(WaveAngle)*sin(WaveAngle))*(gamma*Mach1*Mach1*sin(WaveAngle)*sin(WaveAngle)+1);
    Mach2Star = sqrt(VelocityR1 * VelocityR1 + VelocityPsi1*VelocityPsi1);
    Mach2 = sqrt(2.0/(gamma+1)*Mach2Star*Mach2Star/(1-(gamma-1)/(gamma+1)*Mach2Star*Mach2Star));
    StaticPressure = pow((1+(gamma-1)/2.0*Mach2*Mach2)/(1+(gamma-1)/2.0*Mach1*Mach1), gamma/(gamma-1)) * Pressure;

}




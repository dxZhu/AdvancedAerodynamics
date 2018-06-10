/*
 * main.cpp
 *
 *  Created on: Mar 10, 2017
 *      Author: dxzhu
 */
#include "main.h"

using namespace std;

int main()
{
//number
double Mach1;	//freesream mach number
double Nu;
int N;	//discretized number
double Ni;
double theta0;
double thetaN;
double Mach1Star;	// = V1/a*
double dPsi;
double gamma;
//45 degree case using 4th orderscheme(Runge-Kutta):
arma::colvec Vr_45;
arma::colvec Vpsi_45;
double ConeSemiAngle_45;
arma::colvec Mach_45;
arma::colvec Pressure_45;
arma::colvec Temperature_45;
arma::colvec DensityCone_45;
arma::colvec Cp_45;

//Vector
vector<pair<double,double> > Velocity;
pair<double,double> Vel;
arma::colvec WaveAngle;
arma::colvec ConeSemiAngle;
arma::colvec VelocityR1;	//initial Vr2
arma::colvec VelocityPsi1;	//initial Vpsi2
arma::colvec Pressure1;	//initial pressure ratio p2/p1
arma::colvec StaticPressure1;	//static pressure ratio p02/p01
arma::colvec Temperature1;	//temperature ratio T2/T1
arma::colvec Density1;	//density ratio rou2/rou1
arma::colvec MachCone;	//Mach number at cone surface
arma::colvec PressureCone;
arma::colvec TemperatureCone;
arma::colvec DensityCone;
arma::colvec CpCone;
//Global variable and their temperary store vectors
arma::mat VelocityR;
arma::colvec Vr_T;
arma::mat VelocityPsi;
arma::colvec Vpsi_T;
arma::mat Psi;
arma::colvec Psi_T;
arma::mat Pressure;
arma::colvec P_T;
arma::mat Temperature;
arma::colvec T_T;
arma::mat Density;
arma::colvec D_T;
arma::mat Mach;
arma::colvec M_T;
arma::mat Cp;
arma::colvec Cp_T;

N = 14;	//Discretized wave angles number
Ni = 200;
VelocityR1 = arma::zeros<arma::colvec>(N);
VelocityPsi1 = arma::zeros<arma::colvec>(N);
Pressure1 = arma::zeros<arma::colvec>(N);
StaticPressure1 = arma::zeros<arma::colvec>(N);
Temperature1 = arma::zeros<arma::colvec>(N);
Density1 = arma::zeros<arma::colvec>(N);
MachCone = arma::zeros<arma::colvec>(N);
PressureCone = arma::zeros<arma::colvec>(N);
TemperatureCone = arma::zeros<arma::colvec>(N);
DensityCone = arma::zeros<arma::colvec>(N);
CpCone = arma::zeros<arma::colvec>(N);
VelocityR = arma::zeros<arma::mat>(Ni,N);
VelocityPsi = arma::zeros<arma::mat>(Ni,N);
Vr_T = arma::zeros<arma::colvec>(Ni);
Vpsi_T = arma::zeros<arma::colvec>(Ni);
Psi = arma::zeros<arma::mat>(Ni,N);
Psi_T = arma::zeros<arma::colvec>(Ni);
Pressure = arma::zeros<arma::mat>(Ni,N);
P_T = arma::zeros<arma::colvec>(Ni);
Temperature = arma::zeros<arma::mat>(Ni,N);
T_T = arma::zeros<arma::colvec>(Ni);
Density = arma::zeros<arma::mat>(Ni,N);
D_T = arma::zeros<arma::colvec>(Ni);
Mach = arma::zeros<arma::mat>(Ni,N);
M_T = arma::zeros<arma::colvec>(Ni);
Cp = arma::zeros<arma::mat>(Ni,N);
Cp_T = arma::zeros<arma::colvec>(Ni);
//45 degree
Mach_45 = arma::zeros<arma::colvec>(Ni);
Cp_45 = arma::zeros<arma::colvec>(Ni);
Pressure_45 = arma::zeros<arma::colvec>(Ni);
Temperature_45 = arma::zeros<arma::colvec>(Ni);
DensityCone_45 = arma::zeros<arma::colvec>(Ni);
Vr_45 = arma::zeros<arma::colvec>(Ni);
Vpsi_45 = arma::zeros<arma::colvec>(Ni);

gamma = 1.4;
Mach1 = 2;
Nu = sin(1/Mach1);
theta0 = 0.52333;
thetaN = 1.20367;
Mach1Star = sqrt(((gamma+1)/2 * Mach1 * Mach1)/(1 + (gamma-1)/2 * Mach1 * Mach1));	//M1*
WaveAngle = arma::linspace<arma::colvec>(theta0,thetaN,N);
ConeSemiAngle = arma::zeros<arma::colvec>(N);
cout << "Mach1* is " << Mach1Star << endl;


for(int i=0; i<N; i++){

    InitialCondition( WaveAngle(i),  Mach1,  Mach1Star, Nu, Ni, gamma,
	dPsi, VelocityR1(i), VelocityPsi1(i), Temperature1(i), Density1(i), Pressure1(i), StaticPressure1(i));
//cout << WaveAngle(i) << ' '<< VelocityR1 << ' ' << VelocityPsi1 << endl;


   Iteration(VelocityR1(i), VelocityPsi1(i), gamma, dPsi, Ni, WaveAngle(i),
	   Vr_T, Vpsi_T, Velocity, Psi_T, ConeSemiAngle(i), Vr_45, Vpsi_45, ConeSemiAngle_45);
   VelocityR.col(i) = Vr_T;
   VelocityPsi.col(i) = Vpsi_T;
   Psi.col(i) = Psi_T;


   //Get properties for every point: mach, pressure, temperature, density, Cp
       for(int j=0; j<Ni; j++){
   	Properties(Mach1, gamma, Vr_T(j), Vpsi_T(j), StaticPressure1(i),
   		    M_T(j), P_T(j), T_T(j), D_T(j), Cp_T(j)); }

       cout << "run " << endl;
       Mach.col(i) = M_T;
       Pressure.col(i) = P_T;
       Temperature.col(i) = T_T;
       Density.col(i) = D_T;
       Cp.col(i) = Cp_T;
//Get cone surface properties

cout << "run " << endl;
}

//Get properties for 45degree using R_K:
for(int j=0; j<Ni; j++){
    Properties(Mach1, gamma, Vr_45(j), Vpsi_45(j), StaticPressure1(5),
    	    Mach_45(j), Pressure_45(j), Temperature_45(j), DensityCone_45(j), Cp_45(j));}



//Print solutions
    PrintSolution(WaveAngle, ConeSemiAngle, MachCone, CpCone,
	    ConeSemiAngle_45, Mach_45, Cp_45, Pressure_45, Temperature_45, DensityCone_45, Vr_45, Vpsi_45,
	    VelocityR, VelocityPsi, Psi, Mach, Pressure, Temperature, Density, Cp);

}

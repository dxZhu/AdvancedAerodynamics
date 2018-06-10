/*
 * AxisSymmetry.cpp
 *
 *  Created on: Mar 31, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

void AxisSymmetry(double u1, double v1, double x1, double y1, double gamma, double del, double ep, int Ni, double& x3, double& y3, double& u3, double& v3,
	arma::colvec& Ym, arma::colvec& Yp, arma::colvec& Xm, arma::colvec& Um, arma::colvec& Up, arma::colvec& Vm, arma::colvec& Vp,
	arma::colvec& MachStarm, arma::colvec& Machm,
	arma::colvec& alpham, arma::colvec& thetam, arma::colvec& lamdam,
	arma::colvec& Qm, arma::colvec& Rm, arma::colvec& Sm,
	arma::colvec& X3, arma::colvec& Y3, arma::colvec& U3, arma::colvec& V3)
{
double a;
double dx;
arma::colvec Vel;

Vel = arma::zeros<arma::colvec>(Ni);
Um = arma::zeros<arma::colvec>(Ni);
Up = arma::zeros<arma::colvec>(Ni);
Vm = arma::zeros<arma::colvec>(Ni);
Vp = arma::zeros<arma::colvec>(Ni);
Ym = arma::zeros<arma::colvec>(Ni);
Yp = arma::zeros<arma::colvec>(Ni);
Xm = arma::zeros<arma::colvec>(Ni);
MachStarm = arma::zeros<arma::colvec>(Ni);
Machm = arma::zeros<arma::colvec>(Ni);
alpham = arma::zeros<arma::colvec>(Ni);
thetam = arma::zeros<arma::colvec>(Ni);
lamdam = arma::zeros<arma::colvec>(Ni);
Qm = arma::zeros<arma::colvec>(Ni);
Rm = arma::zeros<arma::colvec>(Ni);
Sm = arma::zeros<arma::colvec>(Ni);
X3 = arma::zeros<arma::colvec>(Ni);
Y3 = arma::zeros<arma::colvec>(Ni);
U3 = arma::zeros<arma::colvec>(Ni);
V3 = arma::zeros<arma::colvec>(Ni);
Um(0) = u1;
Vm(0) = v1;
Vel(0) = sqrt(u1*u1 + v1*v1);
Ym(0) = y1;
dx = 0.0005;

for(int i=0; i<Ni; i++){
a = (gamma+1)/2 - (gamma-1)/2*(Um(i)*Um(i)+Vm(i)*Vm(i));
MachStarm(i) = sqrt(Um(i)*Um(i)+Vm(i)*Vm(i));
Machm(i) = sqrt(2.0/(gamma+1)*MachStarm(i)*MachStarm(i)/(1-(gamma-1)/(gamma+1)*MachStarm(i)*MachStarm(i)));
alpham(i) = asin(1/Machm(i));
thetam(i) = atan(Vm(i)/Um(i));
lamdam(i) = tan(thetam(i) - alpham(i));
Qm(i) = Um(i)*Um(i) - a;
Rm(i) = 2*Um(i)*Vm(i) - Qm(i)*lamdam(i);
Sm(i) = del*a*Vm(i)/Ym(i);
X3(i) = x1 + (i+1) * dx;
Y3(i) = y1 + lamdam(i) * (i+1) * dx;
U3(i) = u1 + (Sm(i)* (i+1) * dx + Rm(i)*v1)/Qm(i);
V3(i) = sqrt(u1*u1 + v1*v1) * sin(thetam(i));

if((Y3(i) < ep) && (V3(i)<ep)){break;}
x3 = X3(i);
y3 = Y3(i);
u3 = U3(i);
v3 = V3(i);

//Corrector
Um(i+1) = (u1 + U3(i))/2;
Vm(i+1) = v1/2;
Ym(i+1) = (y1 + Y3(i))/2;

}
}



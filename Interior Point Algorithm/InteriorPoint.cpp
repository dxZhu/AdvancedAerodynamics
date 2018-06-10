/*
 * InteriorPoint.cpp
 *
 *  Created on: Mar 30, 2017
 *      Author: dxzhu
 */
#include "main.h"

using namespace std;

void InteriorPoint(double u1, double u2, double v1, double v2, double x1, double x2, double y1, double y2, double gamma, double del, double ep, int Ni, double& x3, double& y3, double& u3, double& v3,
	arma::colvec& Ym, arma::colvec& Yp, arma::colvec& Xm, arma::colvec& Xp, arma::colvec& Um, arma::colvec& Up, arma::colvec& Vm, arma::colvec& Vp,
	arma::colvec& MachStarm, arma::colvec& MachStarp, arma::colvec& Machm, arma::colvec& Machp,
	arma::colvec& alpham, arma::colvec& alphap, arma::colvec& thetam, arma::colvec& thetap, arma::colvec& lamdam, arma::colvec& lamdap,
	arma::colvec& Qm, arma::colvec& Qp, arma::colvec& Rm, arma::colvec& Rp, arma::colvec& Sm, arma::colvec& Sp,
	arma::colvec& X3, arma::colvec& Y3, arma::colvec& U3, arma::colvec& V3)
{

Um = arma::zeros<arma::colvec>(Ni);
Up = arma::zeros<arma::colvec>(Ni);
Vm = arma::zeros<arma::colvec>(Ni);
Vp = arma::zeros<arma::colvec>(Ni);
Ym = arma::zeros<arma::colvec>(Ni);
Yp = arma::zeros<arma::colvec>(Ni);
Xm = arma::zeros<arma::colvec>(Ni);
Xp = arma::zeros<arma::colvec>(Ni);
MachStarm = arma::zeros<arma::colvec>(Ni);
MachStarp = arma::zeros<arma::colvec>(Ni);
Machm = arma::zeros<arma::colvec>(Ni);
Machp = arma::zeros<arma::colvec>(Ni);
alpham = arma::zeros<arma::colvec>(Ni);
alphap = arma::zeros<arma::colvec>(Ni);
thetam = arma::zeros<arma::colvec>(Ni);
thetap = arma::zeros<arma::colvec>(Ni);
lamdam = arma::zeros<arma::colvec>(Ni);
lamdap = arma::zeros<arma::colvec>(Ni);
Qm = arma::zeros<arma::colvec>(Ni);
Qp = arma::zeros<arma::colvec>(Ni);
Rm = arma::zeros<arma::colvec>(Ni);
Rp = arma::zeros<arma::colvec>(Ni);
Sm = arma::zeros<arma::colvec>(Ni);
Sp = arma::zeros<arma::colvec>(Ni);
X3 = arma::zeros<arma::colvec>(Ni);
Y3 = arma::zeros<arma::colvec>(Ni);
U3 = arma::zeros<arma::colvec>(Ni);
V3 = arma::zeros<arma::colvec>(Ni);

Um(0) = u1;
Up(0) = u2;
Vm(0) = v1;
Vp(0) = v2;
Ym(0) = y1;
Yp(0) = y2+0.001;

double am;
double ap;

for(int i=0; i< Ni-1; i++){

MachStarm(i) = sqrt(Um(i)*Um(i) + Vm(i)*Vm(i));
MachStarp(i) = sqrt(Up(i)*Up(i) + Vp(i)*Vp(i));
Machm(i) = sqrt(2.0*MachStarm(i)*MachStarm(i)/(gamma+1) / (1.0- ((gamma-1)/(gamma+1)*MachStarm(i)*MachStarm(i))));
Machp(i) = sqrt(2.0*MachStarp(i)*MachStarp(i)/(gamma+1) / (1.0- ((gamma-1)/(gamma+1)*MachStarp(i)*MachStarp(i))));
alpham(i) = asin(1.0/Machm(i));
alphap(i) = asin(1.0/Machp(i));
thetam(i) = atan(Vm(i)/Um(i));
thetap(i) = atan(Vp(i)/Up(i));
lamdam(i) = tan(thetam(i) - alpham(i));
lamdap(i) = tan(thetap(i) + alphap(i));
am = (gamma+1)/2.0 - (gamma-1)*(Um(i)*Um(i)+Vm(i)*Vm(i))/2.0;
ap = (gamma+1)/2.0 - (gamma-1)*(Up(i)*Up(i)+Vp(i)*Vp(i))/2.0;
Qm(i) = Um(i)*Um(i) - am;
Qp(i) = Up(i)*Up(i) - ap;
Rm(i) = 2.0*Um(i)*Vm(i) - Qm(i)*lamdam(i);
Rp(i) = 2.0*Up(i)*Vp(i) - Qp(i)*lamdap(i);
Sm(i) = 1.0*del*am*Vm(i)/Ym(i);
Sp(i) = 1.0*del*ap*Vp(i)/Yp(i);

//use CHAR to get X3 and Y3:
X3(i) = -lamdam(i)/(lamdap(i)-lamdam(i))*x1 + (y1-y2)/(lamdap(i)-lamdam(i))+lamdap(i)/(lamdap(i)-lamdam(i))*x2;
Y3(i) = lamdam(i)*(X3(i)-x1) + y1;

//use COMPA to get U3 and V3
V3(i) = (Qm(i)*Qp(i)*(u2-u1) + Qm(i)*Rp(i)*v2-Qp(i)*Rm(i)*v1 + Qm(i)*Sp(i)*(X3(i)-x2)-Qp(i)*Sm(i)*(X3(i)-x1))/(Qm(i)*Rp(i)-Qp(i)*Rm(i));
U3(i) = u1+Rm(i)/Qm(i)*(v1-V3(i)) + Sm(i)/Qm(i)*(X3(i)-x1);

x3 = X3(i);
y3 = Y3(i);
u3 = U3(i);
v3 = V3(i);

//See whether need further iteration?
if(i > 0){
    if(((abs((V3(i)-V3(i-1))/V3(i)) < ep) && (abs((U3(i)-U3(i-1))/U3(i)) < ep) && (abs((Y3(i)-Y3(i-1))/Y3(i)) < ep) && (abs((X3(i)-X3(i-1))/X3(i)) < ep)))
    {break;}
}

//Corrector
Um(i+1) = (u1 + U3(i))/2;
Up(i+1) = (u2 + U3(i))/2;
Vm(i+1) = (v1 + V3(i))/2;
Vp(i+1) = (v2 + V3(i))/2;
Ym(i+1) = (y1 + Y3(i))/2;
Yp(i+1) = (y2 + Y3(i))/2;

}

}



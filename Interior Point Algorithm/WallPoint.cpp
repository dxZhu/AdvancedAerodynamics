/*
 * WallPoint.cpp
 *
 *  Created on: Apr 27, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

void WallPoint(double u2, double v2, double x2, double y2, double gamma, double del, double ep, double Ni, double& x3, double& y3, double& u3, double& v3)
{
    arma::colvec Xp;
    arma::colvec Yp;
    arma::colvec Up;
    arma::colvec Vp;
    arma::colvec MachStarp;
    arma::colvec Machp;
    arma::colvec alphap;
    arma::colvec thetap;
    arma::colvec lamdap;
    arma::colvec ap;
    arma::colvec Qp;
    arma::colvec Rp;
    arma::colvec Sp;
    arma::colvec TX3;
    arma::colvec TY3;
    arma::colvec TU3;
    arma::colvec TV3;

    //Initialize Value for every loop
    Xp = arma::zeros<arma::colvec>(Ni);
    Yp = arma::zeros<arma::colvec>(Ni);
    Up = arma::zeros<arma::colvec>(Ni);
    Vp = arma::zeros<arma::colvec>(Ni);
    MachStarp = arma::zeros<arma::colvec>(Ni);
    Machp = arma::zeros<arma::colvec>(Ni);
    alphap = arma::zeros<arma::colvec>(Ni);
    thetap = arma::zeros<arma::colvec>(Ni);
    lamdap = arma::zeros<arma::colvec>(Ni);
    ap = arma::zeros<arma::colvec>(Ni);
    Qp = arma::zeros<arma::colvec>(Ni);
    Rp = arma::zeros<arma::colvec>(Ni);
    Sp = arma::zeros<arma::colvec>(Ni);
    TX3 = arma::zeros<arma::colvec>(Ni);
    TY3 = arma::zeros<arma::colvec>(Ni);
    TU3 = arma::zeros<arma::colvec>(Ni);
    TV3 = arma::zeros<arma::colvec>(Ni);

    //set first iteration values
    Xp(0) = x2;
    Yp(0) = y2 + 0.0001;
    Up(0) = u2;
    Vp(0) = v2;

    for(int i=0; i<Ni; i++)
    {
	//Set initial value for C+
	MachStarp(i) = sqrt(Up(i)*Up(i) + Vp(i)*Vp(i));
	Machp(i) = sqrt(2.0*MachStarp(i)*MachStarp(i)/(gamma+1) / (1.0- ((gamma-1)/(gamma+1)*MachStarp(i)*MachStarp(i))));
	alphap(i) = asin(1.0/Machp(i));
	thetap(i) = atan(Vp(i)/Up(i));
	lamdap(i) = tan(thetap(i) + alphap(i));
	ap(i) = (gamma+1)/2.0 - (gamma-1)*(Up(i)*Up(i)+Vp(i)*Vp(i))/2.0;
	Qp(i) = Up(i)*Up(i) - ap(i);
	Rp(i) = 2.0*Up(i)*Vp(i) - Qp(i)*lamdap(i);
	Sp(i) = del*ap(i)*Vp(i)/Yp(i);

	//Use wall function and CHAR to get x3 and y3
	TX3(i) = 1.0/5.0 * (6.0 + 4*Xp(i) - 2*Yp(i) - sqrt(11 - 4.0*Xp(i)*Xp(i) + 4.0*(Yp(i)-3)*Xp(i) + 6*Yp(i) - Yp(i)*Yp(i)));
	TY3(i) = 3.0 - sqrt(4.0 - TX3(i)*TX3(i));

	//Use wall function and COMPA to get u3 and v3
	TU3(i) = (Qp(i)*Up(i) + Rp(i)*Vp(i) + Sp(i)*(TX3(i)-Xp(i)))/(Qp(i) + Rp(i)*TX3(i)/sqrt(4.0 - TX3(i)*TX3(i)));
	TV3(i) = TX3(i)*TU3(i)/sqrt(4.0 - TX3(i)*TX3(i));

	x3 = TX3(i);
	y3 = TY3(i);
	u3 = TU3(i);
	v3 = TV3(i);
/*	cout << "x3 is " << x3;
	cout << "y3 is " << y3;
	cout << "u3 is " << u3;
	cout << "v3 is " << v3 << endl;*/

	//See whether need further iteration?
	if(i > 0){
	    if(((abs((TV3(i)-TV3(i-1))/TV3(i)) < ep) && (abs((TU3(i)-TU3(i-1))/TU3(i)) < ep) && (abs((TY3(i)-TY3(i-1))/TY3(i)) < ep) && (abs((TX3(i)-TX3(i-1))/TX3(i)) < ep)))
	    {break;}
	}

	//Corrector
	Up(i+1) = (u2 + TU3(i))/2;
	Vp(i+1) = (v2 + TV3(i))/2;
	Yp(i+1) = (y2 + TY3(i))/2;
	Xp(i+1) = (x2 + TX3(i))/2;

    }




}



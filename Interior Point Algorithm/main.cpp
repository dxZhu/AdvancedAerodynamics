/*
 * main.cpp
 *
 *  Created on: Mar 30, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

int main()
{
//numbers:
double gamma;
double ep;
double del;
int Ni; //reserve space
int Np;
int Nf;
double x3;
double y3;
double u3;
double v3;
double u1;
double v1;
double x1;
double y1;


Ni = 20000;
Np = 21;
gamma = 1.4;
del = 1;
ep = 0.0001;
Nf = 2*Np -1;

//Initial matrix and vectors:
arma::colvec U1;
arma::colvec V1;
arma::colvec X1;
arma::colvec Y1;
arma::colvec Uac;
arma::colvec Vac;
arma::colvec Xac;
arma::colvec Yac;
arma::colvec Ulm;
arma::colvec Vlm;
arma::colvec Xlm;
arma::colvec Ylm;

//variables
arma::mat thetap;
arma::colvec Tthetap;
arma::mat thetam;
arma::colvec Tthetam;
arma::mat alphap;
arma::colvec Talphap;
arma::mat alpham;
arma::colvec Talpham;
arma::mat lamdap;
arma::colvec Tlamdap;
arma::mat lamdam;
arma::colvec Tlamdam;
arma::mat Up;
arma::colvec TUp;
arma::mat Um;
arma::colvec TUm;
arma::mat Vp;
arma::colvec TVp;
arma::mat Vm;
arma::colvec TVm;
arma::mat Qp;
arma::colvec TQp;
arma::mat Qm;
arma::colvec TQm;
arma::mat Rp;
arma::colvec TRp;
arma::mat Rm;
arma::colvec TRm;
arma::mat Sp;
arma::colvec TSp;
arma::mat Sm;
arma::colvec TSm;
arma::mat Yp;
arma::colvec TYp;
arma::mat Ym;
arma::colvec TYm;
arma::mat Xm;
arma::colvec TXm;
arma::mat Xp;
arma::colvec TXp;
arma::mat Machp;
arma::colvec TMachp;
arma::mat Machm;
arma::colvec TMachm;
arma::mat MachStarp;
arma::colvec TMachStarp;
arma::mat MachStarm;
arma::colvec TMachStarm;
arma::mat X3;
arma::colvec TX3;
arma::mat Y3;
arma::colvec TY3;
arma::mat U3;
arma::colvec TU3;
arma::mat V3;
arma::colvec TV3;
arma::mat Xf;
arma::mat Yf;
arma::mat Uf;
arma::mat Vf;


//Those variables are used to design cacellation domain design:
//Cancellation data
double vpup;
arma::mat Xc;
arma::mat Yc;
arma::mat Uc;
arma::mat Vc;
arma::colvec theta;
theta = arma::zeros<arma::colvec>(Nf-1);
//Cancleation data
Xc = arma::zeros<arma::mat>(Nf,Nf);
Yc = arma::zeros<arma::mat>(Nf,Nf);
Uc = arma::zeros<arma::mat>(Nf,Nf);
Vc = arma::zeros<arma::mat>(Nf,Nf);



string filename1 = "input.dat";
string filename2 = "input_AC.dat";
string filename3 = "input_LM.dat";
string filename4 = "input_CancelLine.dat";
//Calculate ACline:
 x1 = 0.25;
 y1 = 1.0157;
 u1 = 1.2878;
 v1 = 0.15549;

//Initiate all the matrix and vectors will be initiated in interiorPoint

/*Um = arma::zeros<arma::mat>(Ni,Np-1);
Up = arma::zeros<arma::mat>(Ni,Np-1);
Vm = arma::zeros<arma::mat>(Ni,Np-1);
Vp = arma::zeros<arma::mat>(Ni,Np-1);
Ym = arma::zeros<arma::mat>(Ni,Np-1);
Yp = arma::zeros<arma::mat>(Ni,Np-1);
Xm = arma::zeros<arma::mat>(Ni,Np-1);
Xp = arma::zeros<arma::mat>(Ni,Np-1);
MachStarm = arma::zeros<arma::mat>(Ni,Np-1);
MachStarp = arma::zeros<arma::mat>(Ni,Np-1);
Machm = arma::zeros<arma::mat>(Ni,Np-1);
Machp = arma::zeros<arma::mat>(Ni,Np-1);
alpham = arma::zeros<arma::mat>(Ni,Np-1);
alphap = arma::zeros<arma::mat>(Ni,Np-1);
thetam = arma::zeros<arma::mat>(Ni,Np-1);
thetap = arma::zeros<arma::mat>(Ni,Np-1);
lamdam = arma::zeros<arma::mat>(Ni,Np-1);
lamdap = arma::zeros<arma::mat>(Ni,Np-1);
Qm = arma::zeros<arma::mat>(Ni,Np-1);
Qp = arma::zeros<arma::mat>(Ni,Np-1);
Rm = arma::zeros<arma::mat>(Ni,Np-1);
Rp = arma::zeros<arma::mat>(Ni,Np-1);
Sm = arma::zeros<arma::mat>(Ni,Np-1);
Sp = arma::zeros<arma::mat>(Ni,Np-1);*/
X3 = arma::zeros<arma::mat>(Np,2*Np);
Y3 = arma::zeros<arma::mat>(Np,2*Np);
U3 = arma::zeros<arma::mat>(Np,2*Np);
V3 = arma::zeros<arma::mat>(Np,2*Np);
//project forward for design mach number along circle wall
Xf = arma::zeros<arma::mat>(Nf,2*Np);
Yf = arma::zeros<arma::mat>(Nf,2*Np);
Uf = arma::zeros<arma::mat>(Nf,2*Np);
Vf = arma::zeros<arma::mat>(Nf,2*Np);


//ReadInput(filename2, Nf, Uac, Vac, Xac, Yac);
//ReadInput(filename1, Np, U1, V1, X1, Y1);
//ReadInput(filename3, Nf, Ulm, Vlm, Xlm, Ylm);
ReadInput(filename4, Nf, Ulm, Vlm, Xlm, Ylm);

/*//calculate the Domain of dependence !
    //Initial value domain
    X3.col(0) = X1;
    Y3.col(0) = Y1;
    U3.col(0) = U1;
    V3.col(0) = V1;

for(int i=0; i<2*Np-1; i++){

    int m = floor(i/2);
    if(i%2 != 0)
    {
	for(int p=0; p<Np-2-m; p++)
	{

	    InteriorPoint(U3(p+1,i), U3(p,i), V3(p+1,i), V3(p,i), X3(p+1,i), X3(p,i), Y3(p+1,i), Y3(p,i), gamma, del, ep, Ni, x3, y3, u3, v3,
	    		TYm, TYp, TXm, TXp, TUm, TUp, TVm, TVp,
	    		TMachStarm, TMachStarp, TMachm, TMachp,
	    		Talpham, Talphap, Tthetam, Tthetap, Tlamdam, Tlamdap,
	    		TQm, TQp, TRm, TRp, TSm, TSp,
	    		TX3, TY3, TU3, TV3);
	    X3(p+1,i+1) = x3;
	    Y3(p+1,i+1) = y3;
	    U3(p+1,i+1) = u3;
	    V3(p+1,i+1) = v3;
	    cout << x3 << y3 << u3 << v3 << endl;

	}
	AxisSymmetry(U3(0,i), V3(0,i), X3(0,i), Y3(0,i), gamma, del, ep, Ni, x3, y3, u3, v3,
			TYm, TYp, TXm, TUm, TUp, TVm, TVp,
			TMachStarm, TMachm,
			Talpham, Tthetam, Tlamdam,
			TQm, TRm, TSm,
			TX3, TY3, TU3, TV3);
	X3(0,i+1) = x3;
	Y3(0,i+1) = y3;
	U3(0,i+1) = u3;
	V3(0,i+1) = v3;
	cout << x3 << y3 << u3 << v3 << endl;
	cout << "Finish " << i << endl;

    }else{
	for(int p=0; p<Np-m-1; p++)
	{
	    InteriorPoint(U3(p+1,i), U3(p,i), V3(p+1,i), V3(p,i), X3(p+1,i), X3(p,i), Y3(p+1,i), Y3(p,i), gamma, del, ep, Ni, x3, y3, u3, v3,
		    TYm, TYp, TXm, TXp, TUm, TUp, TVm, TVp,
		    TMachStarm, TMachStarp, TMachm, TMachp,
		    Talpham, Talphap, Tthetam, Tthetap, Tlamdam, Tlamdap,
		    TQm, TQp, TRm, TRp, TSm, TSp,
		    TX3, TY3, TU3, TV3);
	    X3(p,i+1) = x3;
	    Y3(p,i+1) = y3;
	    U3(p,i+1) = u3;
	    V3(p,i+1) = v3;
		}
	cout << x3 << y3 << u3 << v3 << endl;
	cout << "Finish " << i << endl;

    }
cout << X3 << endl;
}*/

/*//In line AC data.
for(int j=0; j<2*Np-1; j++){

    int m = floor(j/2);
    if(j%2 != 0){
	X3(Np-1,j) = X3(Np-m-2,j);
	Y3(Np-1,j) = Y3(Np-m-2,j);
	U3(Np-1,j) = U3(Np-m-2,j);
	V3(Np-1,j) = V3(Np-m-2,j);}else{
	X3(Np-1,j) = X3(Np-1-m,j);
	Y3(Np-1,j) = Y3(Np-1-m,j);
	U3(Np-1,j) = U3(Np-1-m,j);
	V3(Np-1,j) = V3(Np-1-m,j);
    }
}*/

/*
//Nozzle wall design y = f(x) = 1 + (x^2)/2, Md = 3.0
//Use this part to iterate get the first two lines
//Use I.V. as x2 y2 u2 v2 to obtain C+ characteristic line
    Xf.col(0) = Xac;
    Yf.col(0) = Yac;
    Uf.col(0) = Uac;
    Vf.col(0) = Vac;
double Md = 3;

for(int i=1; i<7; i++){
    WallPoint(Uf(1,i-1), Vf(1,i-1), Xf(1,i-1), Yf(1,i-1), gamma, del, ep, Ni, x3, y3, u3, v3);
    Xf(0,i) = x3;
    Yf(0,i) = y3;
    Uf(0,i) = u3;
    Vf(0,i) = v3;

    for(int j=0; j<Nf-2; j++){
    InteriorPoint(Uf(j,i), Uf(j+2,i-1), Vf(j,i), Vf(j+2,i-1), Xf(j,i), Xf(j+2,i-1), Yf(j,i), Yf(j+2,i-1), gamma, del, ep, Ni, x3, y3, u3, v3,
	    TYm, TYp, TXm, TXp, TUm, TUp, TVm, TVp,
	    TMachStarm, TMachStarp, TMachm, TMachp,
	    Talpham, Talphap, Tthetam, Tthetap, Tlamdam, Tlamdap,
	    TQm, TQp, TRm, TRp, TSm, TSp,
	    TX3, TY3, TU3, TV3);
    Xf(j+1,i) = x3;
    Yf(j+1,i) = y3;
    Uf(j+1,i) = u3;
    Vf(j+1,i) = v3;
}
cout << "reach" << endl;
    AxisSymmetry(Uf(Nf-2,i), Vf(Nf-2,i), Xf(Nf-2,i), Yf(Nf-2,i), gamma, del, ep, Ni, x3, y3, u3, v3,
	    TYm, TYp, TXm, TUm, TUp, TVm, TVp,
	    TMachStarm, TMachm,
	    Talpham, Tthetam, Tlamdam,
	    TQm, TRm, TSm,
	    TX3, TY3, TU3, TV3);
    Xf(Nf-1,i) = x3;
    Yf(Nf-1,i) = y3;
    Uf(Nf-1,i) = u3;
    Vf(Nf-1,i) = v3;
    cout << Uf << endl;
}
*/
/*//Get Cancellation domain After LM line just two lines
//Initialization got from LM line:
theta(0) = 0.23927;
Xc.col(0) = Xlm;
Yc.col(0) = Ylm;
Uc.col(0) = Ulm;
Vc.col(0) = Vlm;
for(int i=0; i<2; i++){

    for(int it=0; it<Ni; it++){
    //get wall point
    SlopeWallPoint(Uc(1,i), Vc(1,i), Xc(1,i), Yc(1,i), theta(i), gamma, del, ep, Ni, Xc(0,i), Yc(0,i),
	x3, y3, u3, v3);
    Xc(0,i+1) = x3;
    Yc(0,i+1) = y3;
    Uc(0,i+1) = u3;
    Vc(0,i+1) = v3;

for(int j=0; j<Nf-2-i; j++){
cout << "reach1" << endl;
    InteriorPoint(Uc(j,i+1), Uc(j+2,i), Vc(j,i+1), Vc(j+2,i), Xc(j,i+1), Xc(j+2,i), Yc(j,i+1), Yc(j+2,i), gamma, del, ep, Ni, x3, y3, u3, v3,
	    TYm, TYp, TXm, TXp, TUm, TUp, TVm, TVp,
	    TMachStarm, TMachStarp, TMachm, TMachp,
	    Talpham, Talphap, Tthetam, Tthetap, Tlamdam, Tlamdap,
	    TQm, TQp, TRm, TRp, TSm, TSp,
	    TX3, TY3, TU3, TV3);
cout << "reach2" << endl;
Xc(j+1,i+1) = x3;
Yc(j+1,i+1) = y3;
Uc(j+1,i+1) = u3;
Vc(j+1,i+1) = v3;
}

if(Vc(Nf-2-i,i+1) < 0.0002){break;}

//Update theta
theta(i) = theta(i) - atan(Vc(Nf-2-i,i+1)/Uc(Nf-2-i,i+1));
if(theta(i) < 0.0001){theta(i) = 0.00001;}
theta(i+1) = theta(i);
cout << "i is "<< i << endl << "velocity is " << Vc << endl;
    }
}*/

//using finding root technique method to get other lines in cancelation domain:
Xc.col(0) = Xlm;
Yc.col(0) = Ylm;
Uc.col(0) = Ulm;
Vc.col(0) = Vlm;
cout << Xc << endl;
for(int j=0; j<Nf; j++){
for(int i=0; i<37-j; i++){
    FindRoot(Xc(38-j,j),Yc(38-j,j),x3,y3,u3,v3);
    Xc(37-j,j+1) = x3;
    Yc(37-j,j+1) = y3;
    Uc(37-j,j+1) = u3;
    Vc(37-j,j+1) = v3;
cout << "reach" << endl;
    InteriorPoint(Uc(37-j-i,j+1), Uc(37-j-i,j), Vc(37-j-i,j+1), Vc(37-j-i,j), Xc(37-j-i,j+1), Xc(37-j-i,j), Yc(37-j-i,j+1), Yc(37-j-i,j), gamma, del, ep, Ni, x3, y3, u3, v3,
    	    TYm, TYp, TXm, TXp, TUm, TUp, TVm, TVp,
    	    TMachStarm, TMachStarp, TMachm, TMachp,
    	    Talpham, Talphap, Tthetam, Tthetap, Tlamdam, Tlamdap,
    	    TQm, TQp, TRm, TRp, TSm, TSp,
    	    TX3, TY3, TU3, TV3);
        Xc(36-j-i,j+1)= x3;
        Yc(36-j-i,j+1)= y3;
        Uc(36-j-i,j+1)= u3;
        Vc(36-j-i,j+1)= v3;
    }
}


cout << Xc << endl << Yc << endl << Uc << endl << Vc;
//Print solutions:
PrintSolutions( Xc, Yc, Uc, Vc);
}



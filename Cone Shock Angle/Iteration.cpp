/*
 * Iteration.cpp
 *
 *  Created on: Mar 10, 2017
 *      Author: dxzhu
 */
#include "main.h"

using namespace std;

void Iteration(double VelocityR1, double VelocityPsi1, double gamma, double dPsi, double Ni, double WaveAngle,
	arma::colvec& Vr, arma::colvec& Vpsi, vector<pair<double,double> > Velocity, arma::colvec& Psi,
	double& ConeSemiAngle, arma::colvec& Vr_45, arma::colvec& Vpsi_45, double& ConeSemiAngle_45)
{
double k1;
double k2;
double k3;
double k4;
double l1;
double l2;
double l3;
double l4;
double SoundSpeed; //	a/a*^2
double ep1;
pair<double,double> Vel;

Vr = arma::zeros<arma::colvec>(Ni);
Vpsi = arma::zeros<arma::colvec>(Ni);
Psi = arma::zeros<arma::colvec>(Ni);

ep1 = 1.0;
Vr(0) = VelocityR1;	//initial Vr
Vpsi(0) = VelocityPsi1;	//initial Vpsi
SoundSpeed = (gamma+1)/2 - (gamma-1)/2*(VelocityR1*VelocityR1+VelocityPsi1*VelocityPsi1);

for(int i=0; i<Ni; i++){

	int n=0;	//iteration number for each delta Psi
	Psi(i) = WaveAngle + i*dPsi;
	SoundSpeed = (gamma+1)/2 - (gamma-1)/2*(Vr(i)*Vr(i)+Vpsi(i)*Vpsi(i));
	k1 = Vpsi(i);
	l1 = SoundSpeed*(Vr(i)+Vpsi(i)/tan(Psi(i)))/(Vpsi(i)*Vpsi(i)-SoundSpeed)-Vr(i);
	Vr(i+1) = Vr(i) + dPsi*k1;	//starting point of iteration
	Vpsi(i+1) = Vpsi(i) + dPsi*l1;	//starting point of iterarion
	//do the iteration
	do{
	    SoundSpeed = (gamma+1)/2 - (gamma-1)/2*(Vr(i+1)*Vr(i+1)+Vpsi(i+1)*Vpsi(i+1));
	    k2 = Vpsi(i+1);
	    l2 = SoundSpeed*(Vr(i+1)+Vpsi(i+1)/tan(Psi(i)))/(Vpsi(i+1)*Vpsi(i+1)-SoundSpeed)-Vr(i+1);
	    ep1 =Vr(i+1);

	    Vr(i+1) = Vr(i) + dPsi/2 * (k1 + k2);
	    Vpsi(i+1) = Vpsi(i) + dPsi/2 * (l1 + l2);
	    ep1 = abs(Vr(i+1)-ep1)/Vr(i+1);
	    //cout << n << " Vpsi is " << Vpsi(i+1) << " Vr(i+1) and Vr(i) is " << Vr(i+1) << ' ' << Vr(i) << "Cone angle is " << ConeSemiAngle << endl;
	    n++;
		}while(ep1>0.0000001);

	ConeSemiAngle =WaveAngle + i*dPsi;
	if(Vpsi(i+1)>0){Vr(i+1) = ((abs(Vpsi(i))<abs(Vpsi(i+1))) ? Vr(i):Vr(i+1));
			Vpsi(i+1) = ((abs(Vpsi(i))<abs(Vpsi(i+1))) ? Vpsi(i):Vpsi(i+1));
			Vel.first = Vr(i+1);
			Vel.second = Vpsi(i+1);
			break;}

    }
    Velocity.push_back(Vel);

	//Print out for 45 degree
	if(WaveAngle==(1.20367-0.52333)/13*5+0.52333){
	    ep1 = 1.0;
	    Vr_45(0) = VelocityR1;	//initial Vr
	    Vpsi_45(0) = VelocityPsi1;	//initial Vpsi
	    SoundSpeed = (gamma+1)/2 - (gamma-1)/2*(VelocityR1*VelocityR1+VelocityPsi1*VelocityPsi1);
	    //4th order scheme
	    for(int i=0; i<Ni; i++){
	    	int n=0;	//iteration number for each delta Psi
	    	Psi(i) = WaveAngle + i*dPsi;
	    	SoundSpeed = (gamma+1)/2 - (gamma-1)/2*(Vr_45(i)*Vr_45(i)+Vpsi_45(i)*Vpsi_45(i));
	    	k1 = Vpsi_45(i);
	    	l1 = SoundSpeed*(Vr_45(i)+Vpsi_45(i)/tan(Psi(i)))/(Vpsi_45(i)*Vpsi_45(i)-SoundSpeed)-Vr_45(i);
	    	k2 = Vpsi_45(i+1) + l1*dPsi/2;
	    	l2 = SoundSpeed*((Vr_45(i+1)+k1*dPsi/2) + (Vpsi_45(i+1)+l1*dPsi/2)/tan(Psi(i)))/((Vpsi_45(i+1)+l1*dPsi/2)*(Vpsi_45(i+1)+l1*dPsi/2)-SoundSpeed)-(Vr_45(i+1)+k1*dPsi/2);
	    	k3 = Vpsi_45(i+1) + l2*dPsi/2;
	    	l3 = SoundSpeed*((Vr_45(i+1)+k2*dPsi/2) + (Vpsi_45(i+1)+l2*dPsi/2)/tan(Psi(i)))/((Vpsi_45(i+1)+l2*dPsi/2)*(Vpsi_45(i+1)+l2*dPsi/2)-SoundSpeed)-(Vr_45(i+1)+k2*dPsi/2);
	    	k4 = Vpsi_45(i+1) + l3*dPsi/2;;
	    	l4 = SoundSpeed*((Vr_45(i+1)+k3*dPsi/2) + (Vpsi_45(i+1)+l3*dPsi/2)/tan(Psi(i)))/((Vpsi_45(i+1)+l3*dPsi/2)*(Vpsi_45(i+1)+l3*dPsi/2)-SoundSpeed)-(Vr_45(i+1)+k3*dPsi/2);
	    	Vr_45(i+1) = Vr_45(i) + dPsi/6*(k1+2*k2+2*k3+k4);	//starting point of iteration
	    	Vpsi_45(i+1) = Vpsi_45(i) + dPsi/6*(l1+2*l2+2*l3+l4);	//starting point of iterarion
	    	//do the iteration
	    	do{
	    	SoundSpeed = (gamma+1)/2 - (gamma-1)/2*(Vr_45(i+1)*Vr_45(i+1)+Vpsi_45(i+1)*Vpsi_45(i+1));
	    	k2 = Vpsi_45(i+1);
	    	l2 = SoundSpeed*(Vr_45(i+1)+Vpsi_45(i+1)/tan(Psi(i)))/(Vpsi_45(i+1)*Vpsi_45(i+1)-SoundSpeed)-Vr_45(i+1);
	    	ep1 =Vr_45(i+1);

	    	Vr_45(i+1) = Vr_45(i) + dPsi/2 * (k1 + k2);
	    	Vpsi_45(i+1) = Vpsi_45(i) + dPsi/2 * (l1 + l2);
	    	ep1 = abs(Vr_45(i+1)-ep1)/Vr_45(i+1);
	    	    //cout << n << " Vpsi is " << Vpsi(i+1) << " Vr(i+1) and Vr(i) is " << Vr(i+1) << ' ' << Vr(i) << "Cone angle is " << ConeSemiAngle << endl;
	    	    n++;
	    		}while(ep1>0.0000001);

	    	ConeSemiAngle_45 =WaveAngle + i*dPsi;
	    	if(Vpsi_45(i+1)>0){Vr_45(i+1) = ((abs(Vpsi_45(i))<abs(Vpsi_45(i+1))) ? Vr_45(i):Vr_45(i+1));
	    			Vpsi_45(i+1) = ((abs(Vpsi_45(i))<abs(Vpsi_45(i+1))) ? Vpsi_45(i):Vpsi_45(i+1));
	    			break;}

	        }
	}



}




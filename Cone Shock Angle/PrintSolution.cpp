/*
 * PrintSolution.cpp
 *
 *  Created on: Mar 12, 2017
 *      Author: dxzhu
 */
#include "main.h"

using namespace std;

void PrintSolution(arma::colvec& WaveAngle, arma::colvec& ConeSemiAngle,  arma::colvec& MachCone,   arma::colvec& CpCone,
	double& ConeSemiAngle_45, arma::colvec& MachCone_45, arma::colvec& CpCone_45, arma::colvec& Pressure_45, arma::colvec& Temperature_45, arma::colvec& DensityCone_45,arma::colvec& Vr_45, arma::colvec& Vpsi_45,
	arma::mat& VelocityR, arma::mat& VelocityPsi, arma::mat& Psi, arma::mat& Mach, arma::mat& Pressure, arma::mat& Temperature, arma::mat& Density, arma::mat& Cp)
{

  ofstream fout;

  fout.open ("WaveAngle.dat");
  fout << WaveAngle;
  fout.close();

  fout.open ("ConeSemiAngle.dat");
  fout << ConeSemiAngle;
  fout.close();

  fout.open ("MachCone.dat");
  fout << MachCone;
  fout.close();

  fout.open ("CpCone.dat");
  fout << CpCone;
  fout.close();

  fout.open ("PrintOut45_4th.dat");
  fout << "Vr_45 and Vpsi is : " << endl
	  << Vr_45 <<' '<<Vpsi_45<< endl;
  fout.close();

  fout.open ("VelocityR.dat");
  fout << VelocityR;
  fout.close();

  fout.open ("VelocityPsi.dat");
  fout << VelocityPsi;
  fout.close();

  fout.open ("Psi.dat");
  fout << Psi;
  fout.close();

  fout.open ("Mach.dat");
  fout << Mach;
  fout.close();

  fout.open ("Pressure.dat");
  fout << Pressure;
  fout.close();

  fout.open ("Temperature.dat");
  fout << Temperature;
  fout.close();


  fout.open ("Density.dat");
  fout << Density;
  fout.close();

  fout.open ("Cp.dat");
  fout << Cp;
  fout.close();
      /*cout<<"Wave Angle is \n";
  cout<<WaveAngle<<endl;
  cout<<"Mach number at cone surface is \n";
  cout<<MachCone<<endl;
  cout<<"ConeSemi Angle is \n";
  cout<<ConeSemiAngle<<endl;*/

}


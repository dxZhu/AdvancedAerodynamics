/*
 * main.h
 *
 *  Created on: Mar 10, 2017
 *      Author: dxzhu
 */

#ifndef MAIN_H_
#define MAIN_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// INCLUDE ARMADILLO LIBRARY DEFINITIONS
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <math.h>

using namespace std;


void InitialCondition(double WaveAngle, double Mach1, double Mach1Star, double Nu, double Ni, double gamma,
	double& dPsi, double& VelocityR1, double& VelocityPsi1, double& Temperature, double& Density, double& Pressure, double& StaticPressure);


void Iteration(double VelocityR1, double VelocityPsi1, double gamma, double dPsi, double Ni, double WaveAngle,
	arma::colvec& Vr, arma::colvec& Vpsi, vector<pair<double,double> > Velocity, arma::colvec& Psi,
	double& ConeSemiAngle, arma::colvec& Vr_45, arma::colvec& Vpsi_45, double& ConeSemiAngle_45);


void Properties(double Mach1, double gamma, double Vr, double Vpsi, double StaticPressure,
	double& Mach, double& Pressure, double& Temperature, double& Density, double& Cp);


void PrintSolution(arma::colvec& WaveAngle, arma::colvec& ConeSemiAngle,  arma::colvec& MachCone,   arma::colvec& CpCone,
	double& ConeSemiAngle_45, arma::colvec& MachCone_45, arma::colvec& CpCone_45, arma::colvec& Pressure_45, arma::colvec& Temperature_45, arma::colvec& DensityCone_45, arma::colvec& Vr_45, arma::colvec& Vpsi_45,
	arma::mat& VelocityR, arma::mat& VelocityPsi, arma::mat& Psi, arma::mat& Mach, arma::mat& Pressure, arma::mat& Temperature, arma::mat& Density, arma::mat& Cp);

#endif /* MAIN_H_ */

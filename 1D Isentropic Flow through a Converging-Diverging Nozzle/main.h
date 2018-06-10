/*
 * main.h
 *
 *  Created on: May 2, 2017
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

void GetInitialValuePlane(double x0In, double x0Out, double rho0In, double rho0Out, double p0In, double p0Out, double v0In, double v0Out, double T0In, double T0Out,
	double dx, int spacestep, int timestep,
	arma::rowvec& area, arma::rowvec& space, arma::mat& velocity, arma::mat& pressure, arma::mat& density, arma::mat& temperature);

void GetPredictor(double dx, double dt, int spacestep, int n, double gamma,
	arma::rowvec area, arma::mat& velocity, arma::mat& pressure, arma::mat& density, arma::mat& temperature);

void GetCorrector(double dx, double dt, int spacestep, int n, double gamma,
	arma::rowvec area, arma::mat& velocity, arma::mat& pressure, arma::mat& density, arma::mat& temperature);

void PrintSolutions(arma::rowvec area, arma::rowvec space, arma::mat velocity, arma::mat pressure, arma::mat density, arma::mat temperature);

void ExtractSolutions(arma::rowvec area, arma::rowvec space, arma::mat velocity, arma::mat pressure, arma::mat density, arma::mat temperature);

#endif /* MAIN_H_ */

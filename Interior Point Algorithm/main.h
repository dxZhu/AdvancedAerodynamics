/*
 * main.h
 *
 *  Created on: Mar 30, 2017
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

void ReadInput(const std::string& filename, int Np,
	arma::colvec& U1, arma::colvec& V1, arma::colvec& X1, arma::colvec& Y1);

void InteriorPoint(double u1, double u2, double v1, double v2, double x1, double x2, double y1, double y2, double gamma, double del, double ep, int Ni, double& x3, double& y3, double& u3, double& v3,
	arma::colvec& Ym, arma::colvec& Yp, arma::colvec& Xm, arma::colvec& Xp, arma::colvec& Um, arma::colvec& Up, arma::colvec& Vm, arma::colvec& Vp,
	arma::colvec& MachStarm, arma::colvec& MachStarp, arma::colvec& Machm, arma::colvec& Machp,
	arma::colvec& alpham, arma::colvec& alphap, arma::colvec& thetam, arma::colvec& thetap, arma::colvec& lamdam, arma::colvec& lamdap,
	arma::colvec& Qm, arma::colvec& Qp, arma::colvec& Rm, arma::colvec& Rp, arma::colvec& Sm, arma::colvec& Sp,
	arma::colvec& X3, arma::colvec& Y3, arma::colvec& U3, arma::colvec& V3);

void AxisSymmetry(double u1, double v1, double x1, double y1, double gamma, double del, double ep, int Ni, double& x3, double& y3, double& u3, double& v3,
	arma::colvec& Ym, arma::colvec& Yp, arma::colvec& Xm, arma::colvec& Um, arma::colvec& Up, arma::colvec& Vm, arma::colvec& Vp,
	arma::colvec& MachStarm, arma::colvec& Machm,
	arma::colvec& alpham, arma::colvec& thetam, arma::colvec& lamdam,
	arma::colvec& Qm, arma::colvec& Rm, arma::colvec& Sm,
	arma::colvec& X3, arma::colvec& Y3, arma::colvec& U3, arma::colvec& V3);

void WallPoint(double u2, double v2, double x2, double y2, double gamma, double del, double ep, double Ni, double& x3, double& y3, double& u3, double& v3);

void SlopeWallPoint(double u2, double v2, double x2, double y2, double theta0, double gamma, double del, double ep, double Ni, double xw, double yw,
	double& x3, double& y3, double& u3, double& v3);

void FindRoot(double xp, double yp, double& x3, double& y3, double& u3, double& v3);

void PrintSolutions(arma::mat X3, arma::mat Y3, arma::mat U3, arma::mat V3);


#endif /* MAIN_H_ */

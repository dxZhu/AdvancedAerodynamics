/*
 * ExtractSolutions.cpp
 *
 *  Created on: May 4, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

void ExtractSolutions(arma::rowvec area, arma::rowvec space, arma::mat velocity, arma::mat pressure, arma::mat density, arma::mat temperature)
{

    ofstream fout;

    fout.open ("InitialValuePlane.dat");
    fout << "Initial value data for the first step" << endl;
    fout << space;
    fout << area;
    fout << velocity.row(0);
    fout << pressure.row(0);
    fout << density.row(0);
    fout << temperature.row(0);
    fout.close();

    fout.open ("100thValuePlane.dat");
    fout << "Initial value data for the first step" << endl;
    fout << space;
    fout << area;
    fout << velocity.row(99);
    fout << pressure.row(99);
    fout << density.row(99);
    fout << temperature.row(99);
    fout.close();

    fout.open ("200thValuePlane.dat");
    fout << "Initial value data for the first step" << endl;
    fout << space;
    fout << area;
    fout << velocity.row(199);
    fout << pressure.row(199);
    fout << density.row(199);
    fout << temperature.row(199);
    fout.close();

    fout.open ("300thValuePlane.dat");
    fout << space;
    fout << area;
    fout << velocity.row(299);
    fout << pressure.row(299);
    fout << density.row(299);
    fout << temperature.row(299);
    fout.close();

    fout.open ("400thValuePlane.dat");
    fout << space;
    fout << area;
    fout << velocity.row(399);
    fout << pressure.row(399);
    fout << density.row(399);
    fout << temperature.row(399);
    fout.close();

    fout.open ("500thValuePlane.dat");
    fout << space;
    fout << area;
    fout << velocity.row(499);
    fout << pressure.row(499);
    fout << density.row(499);
    fout << temperature.row(499);
    fout.close();

    fout.open ("PressureWithTime.dat");
    for(int i=0; i<500; i++){
    fout << i+1 << '	' << pressure(i,50) << '	' << velocity(i,50) << '	' << density(i,50) << '	' << temperature(i,50) << endl;}
    fout.close();

}

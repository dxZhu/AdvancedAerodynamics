/*
 * PrintSolutions.cpp
 *
 *  Created on: May 2, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

void PrintSolutions(arma::rowvec area, arma::rowvec space, arma::mat velocity, arma::mat pressure, arma::mat density, arma::mat temperature)
{

    ofstream fout;

    fout.open ("area.dat");
    fout << area;
    fout.close();

    fout.open ("space.dat");
    fout << space;
    fout.close();

    fout.open ("velocity.dat");
    fout << velocity;
    fout.close();

    fout.open ("pressure.dat");
    fout << pressure;
    fout.close();

    fout.open ("density.dat");
    fout << density;
    fout.close();

    fout.open ("temperature.dat");
    fout << temperature;
    fout.close();
/*
    fout.open ("Um.dat");
    fout << Um;
    fout.close();

    fout.open ("Up.dat");
    fout << Up;
    fout.close();

    fout.open ("Vm.dat");
    fout << Vm;
    fout.close();

    fout.open ("Vp.dat");
    fout << Vp;
    fout.close();

*/



}





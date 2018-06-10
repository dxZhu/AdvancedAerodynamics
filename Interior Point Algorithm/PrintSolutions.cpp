/*
 * PrintSolutions.cpp
 *
 *  Created on: Mar 30, 2017
 *      Author: dxzhu
 */
#include "main.h"

using namespace std;

void PrintSolutions(arma::mat X3, arma::mat Y3, arma::mat U3, arma::mat V3)
{

    ofstream fout;

    fout.open ("U3.dat");
    fout << U3;
    fout.close();

    fout.open ("V3.dat");
    fout << V3;
    fout.close();

    fout.open ("X3.dat");
    fout << X3;
    fout.close();

    fout.open ("Y3.dat");
    fout << Y3;
    fout.close();
/*
    fout.open ("Ym.dat");
    fout << Ym;
    fout.close();

    fout.open ("Yp.dat");
    fout << Yp;
    fout.close();

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


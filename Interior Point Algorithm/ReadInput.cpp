/*
 * ReadInput.cpp
 *
 *  Created on: Apr 6, 2017
 *      Author: dxzhu
 */

#include "main.h"

using namespace std;

void ReadInput(const std::string& filename, int Np,
	arma::colvec& U1, arma::colvec& V1, arma::colvec& X1, arma::colvec& Y1){

ifstream fout(filename.c_str());

U1 = arma::zeros<arma::colvec>(Np);
V1 = arma::zeros<arma::colvec>(Np);
X1 = arma::zeros<arma::colvec>(Np);
Y1 = arma::zeros<arma::colvec>(Np);

for(int i=0; i<Np; i++){

    fout >> X1(i) >> Y1(i) >> U1(i) >> V1(i);
}

fout.close();

}



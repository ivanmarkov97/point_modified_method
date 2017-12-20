// ConsoleApplication2.cpp: îïðåäåëÿåò òî÷êó âõîäà äëÿ êîíñîëüíîãî ïðèëîæåíèÿ.
//

//#include "stdafx.h"
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
//#include "stdafx.h"

using namespace std;

/*double tp = 0;
double dtp = 1e-5;
double t = 0.0;
int ti = 0;
double dt = 0.0;
double tk = 1e-3;
double R10 = 1000.0;
double E10 = 10.0;
double T10 = 1e-4;
double W10 = 2.0*3.14 / T10;
double E9 = 5.0;
double R12 = 1000.0;
double L6 = 0.00001;
double C7 = 1e-6;
double R9 = 1000.0;
double C6 = 1e-6;
double Rb = 20.0;
double Rd = 1000000.0;
double It = 1e-12;
double MFt = 0.026;
double Cb = 2e-12;
double sigma = 0.01;
double R = 1.0 ;
int dem = 4;*/

double tp = 0.0;
double dtp = 1.0e-5;
double t = 0.0;
int ti = 0;
double dt = 0.0;
double tk = 1.0e-3;
double T10 = 1.0e-3;
double W10 = 2.0*3.14 / T10;
double Ea = 1.0;
double Cb1 = 2.0e-12;
double Cb2 = 2.0e-12;
double C1 = 1.0e-6;
double C2 = 1.0e-6;
double Rb1 = 20.0;
double Rb2 = 20.0;
double Ru1 = 1.0e+6;
double Ru2 = 1.0e+6;
double R = 1.0;
double L = 2.53e-4;
double It = 1.0e-12;
double Mft = 0.026;
int dim = 5;
double sigma = 0.01;


vector<double> count_V_P(vector<double> V_n, vector<double> V_P) {
	double E = Ea*sin(W10 * t);
	V_P[0] =  ((E - V_n[0] + V_n[2])/Rb1 + (E - V_n[0] - V_n[1] - V_n[3])/Rb2 - V_n[4])/C1;
	V_P[1] =  ((E - V_n[0] - V_n[1] - V_n[3])/Rb2 - V_n[1]/R)/C2;
	V_P[2] =  (-(E - V_n[0] + V_n[2])/Rb1 - V_n[2]/Ru1 - It*(exp(V_n[2]/Mft) - 1))/Cb1;
	V_P[3] =  ((E - V_n[0] - V_n[1] - V_n[3])/Rb2 - V_n[3]/Ru2 - It*(exp(V_n[3]/Mft) - 1))/Cb2;
	V_P[4] =  V_n[0]/L;
	return V_P;

}

double Norm(vector<double> V_P) {
    double x = 0.0;
		//int check ;
    for(int i = 0 ; i < V_P.size() ; i++){
			x += V_P[i] * V_P[i] ;
			//std::cout << V_P[i] << " ";
		}
		//std::cout << std::endl;
		//cin >> check;
    x = sqrt(x) ;
    return x ;
}

void mps() {
	ofstream fout("result.txt");
  int check;
	vector<double> V(dim);
	vector<double> V_P(dim);
	vector<double> V_n(dim);
	for (int i = 0; i < dim; i++)
		V_n[i] = 1e-8;
	while (t < tk) {
		V_P = count_V_P(V_n, V_P);
		dt = sigma / Norm(V_P);
		for(int i = 0 ; i < V_P.size() ; i ++)
            V[i] = V_P[i] * dt + V_n[i];
		t = t + dt;
		if (t > tp) {
      fout << V[1] << " " ;
      fout << endl ;
			tp += dtp;
		}
		V_n = V;

	}
  fout.close();
}

int main() {
	mps();
    return 0;
}

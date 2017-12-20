// ConsoleApplication2.cpp: îïðåäåëÿåò òî÷êó âõîäà äëÿ êîíñîëüíîãî ïðèëîæåíèÿ.
//

//#include "stdafx.h"
#include<vector>




#include<stdlib.h>
#include<iostream>
#include<math.h>
#include <fstream>
#include<string>
//#include "stdafx.h"

using namespace std;

double tp = 0;
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
int dem = 4;






vector<double> count_V_P(vector<double> V_n, vector<double> V_P) {
	double E10_s = E10*sin(W10 * t);
	V_P[0] =  ((-V_n[0]) / R10 - (-E10_s + E9 + V_n[0] - V_n[2]) / Rb + V_n[3]) / C6 ;
	V_P[1] =  ( - V_n[1] / R9 + V_n[3]) / C7;
	V_P[2] = (( - E10_s + E9 + V_n[0] - V_n[2]) / Rb - V_n[2] / Rd + It*(exp(-V_n[2]/MFt) - 1)) / Cb ;
	V_P[3] = (E10_s - V_n[0] - V_n[1] - V_n[3] * R12) / L6;
	return V_P;

}

double Norm(vector<double> V_P) {
    double x = 0.0;
    for(int i = 0 ; i < V_P.size() ; i ++)
        x += V_P[i] * V_P[i] ;
    x = sqrt(x) ;
    return x ;
}


// c6 c7 cd il


void mps() {
	ofstream fout("result.txt");
	vector<double> V(dem);
	vector<double> V_P(dem);
	vector<double> V_n(dem);
	for (int i = 0; i < dem; i++)
		V_n[i] = 0.0 ;
	while (t < tk) {
		V_P = count_V_P(V_n, V_P);
		dt = sigma / Norm(V_P) ;
//cout << dt << endl ;
		for(int i = 0 ; i < V_P.size() ; i ++)
            V[i] = V_P[i] * dt + V_n[i];
		t = t + dt;
		if (t > tp) {
			fout << t << " ";
			for(int i = 0 ; i < V.size() ; i ++)
                fout << V[i] << " " ;
            fout << endl ;
			tp += dtp;
		}
		V_n = V;

	}

}

int main()
{
	mps();
    return 0;
}

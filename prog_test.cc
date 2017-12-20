#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#define It1 1e-12
#define It2 1e-12
#define Cb1 2e-12
#define Cb2 2e-12
#define Mft 0.026
#define Ru 1000000.0
#define RB 20.0
#define t_min 1.0e-2
#define C1 1.0e-6
#define C2 1.0e-6
#define R 1.0e+3
#define L 2.53e-4
#define E(t) 1.0 * sin(2 * M_PI * t)
#define A1(F4, F2) 1.0 / Mft * It1 * (exp((F4 - F2) / Mft) )
#define A2(F4, F2) -1.0 / Mft * It1 * (exp((F4 - F2) / Mft) )
#define A3(F5, F3) 1.0 / Mft * It2 * (exp((F5 - F3) / Mft) )
#define A4(F5, F3) -1.0 / Mft * It2 * (exp((F5 - F3) / Mft) )
#define N 16
#define TIME 1.0
#define delta 1.0e-2
#define MAX_ITER_STEPS 100
#define EPS1 0.01
#define EPS2 0.04
#define INTEGR_DELTA 0.01

using namespace std;

double t = 0.0;
double dt_prev = 0.0;
double dt = t_min;
double prevStep[N] = {0.0};
double curApprox[N] = {0.0};
double U_m_1[N] = {0.0};
double U_m_2[N] = {0.0};
double eps[N] = {0.0};

string str[N] = {"_UC1", "_UC2", "_UCb1", "_UCb2", "_Il",
  "UC1", "UC2", "UCb1", "UCb2", "Il",
  "fi1", "fi2", "fi3", "fi4", "fi5",
  "Ie"
};

enum BASE {
  _UC1, _UC2, _UCb1, _UCb2, _Il,
  UC1, UC2, UCb1, UCb2, Il,
  fi1, fi2, fi3, fi4, fi5,
  Ie
};

class Scheme{
public:
  Scheme(): size(0), global_matx(NULL), free_part(NULL){};
  Scheme(int n);
  int get_size(){ return size; };
  double **get_global_matx(){ return global_matx; };
  double *get_free_part(){ return free_part; };
  void init_global_matx();
  void init_free_part();
  void show_global_matx();
  void show_free_part();
  friend double* gausse();
  void complete_iter(double*, double* , double*);

private:
  double **global_matx;
  double *free_part;
  int size;
};

Scheme::Scheme(int n){
  size = n;
  free_part = new double[size];
  global_matx = new double*[size];
  for(int i = 0; i < size; i++)
    global_matx[i] = new double[size];
}

void Scheme::init_global_matx(){
  for(int i = 0; i < 5; i++){
    global_matx[i][i] = 1.0;
  }

  for(int i = 0; i < 5; i++){
    global_matx[i][i + 5] = -1.0 / dt;
  }

  for(int i = 5; i < 9; i++){
    global_matx[i][i] = 1.0;
  }

  global_matx[5][10] = global_matx[9][10] = -1.0;
  global_matx[5][11] = global_matx[9][11] = 1.0;
  global_matx[7][11] = global_matx[8][12] = 1.0;
  global_matx[6][12] = global_matx[7][13] = global_matx[8][14] = -1.0;

  global_matx[9][4] = L;

  global_matx[10][0] = C1;
  global_matx[11][0] = -C1;
  global_matx[11][2] = -Cb1;
  global_matx[12][1] = C2;
  global_matx[12][3] = -Cb2;
  global_matx[13][2] = Cb1;
  global_matx[14][3] = Cb2;

  global_matx[10][9] = 1.0;
  global_matx[11][9] = -1.0;

  global_matx[11][11] = -A2(curApprox[BASE::fi4], curApprox[BASE::fi2]) + 1.0 / Ru + 1.0 / RB;
  global_matx[11][13] = - A1(curApprox[BASE::fi4], curApprox[BASE::fi2]) - 1.0 / Ru;
  global_matx[11][14] = -1.0 / RB;
  global_matx[12][12] = -A4(curApprox[BASE::fi5], curApprox[BASE::fi3]) + 1.0 / Ru + 1.0 / R;
  global_matx[12][14] = -1.0 / Ru - A3(curApprox[BASE::fi5], curApprox[BASE::fi3]);
  global_matx[13][11] = -1.0 / Ru + A2(curApprox[BASE::fi4], curApprox[BASE::fi2]);
  global_matx[13][13] = A1(curApprox[BASE::fi4], curApprox[BASE::fi2]) + 1.0 / RB + 1.0 / Ru;
  global_matx[14][11] = -1.0 / RB;
  global_matx[14][12] = -1.0 / Ru + A4(curApprox[BASE::fi5], curApprox[BASE::fi3]);
  global_matx[14][14] = A3(curApprox[BASE::fi5], curApprox[BASE::fi3]) + 1.0 / RB + 1.0 / Ru;

  global_matx[10][15] = 1.0;

  global_matx[15][10] = 1.0;
}

void Scheme::init_free_part(){
  free_part[BASE::_UC1] = - (curApprox[BASE::_UC1] - (curApprox[BASE::UC1] - prevStep[BASE::UC1]) / dt);
  free_part[BASE::_UC2] = - (curApprox[BASE::_UC2] - (curApprox[BASE::UC2] - prevStep[BASE::UC2]) / dt);
  free_part[BASE::_UCb1] = - (curApprox[BASE::_UCb1] - (curApprox[BASE::UCb1] - prevStep[BASE::UCb1]) / dt);
  free_part[BASE::_UCb2] = - (curApprox[BASE::_UCb2] - (curApprox[BASE::UCb2] - prevStep[BASE::UCb2]) / dt);
  free_part[BASE::_Il] = - (curApprox[BASE::_Il] - (curApprox[BASE::Il] - prevStep[BASE::Il]) / dt);
  free_part[BASE::UC1] = - (curApprox[BASE::UC1] - (curApprox[BASE::fi1] - curApprox[BASE::fi2]));
  free_part[BASE::UC2] = - (curApprox[BASE::UC2] - curApprox[BASE::fi3]);
  free_part[BASE::UCb1] = - (curApprox[BASE::UCb1] - (curApprox[BASE::fi4] - curApprox[BASE::fi2]));
  free_part[BASE::UCb2] = - (curApprox[BASE::UCb2] - (curApprox[BASE::fi5] - curApprox[BASE::fi3]));
  free_part[BASE::Il] = - (L * curApprox[BASE::_Il] - (curApprox[BASE::fi1] - curApprox[BASE::fi2]));
  free_part[BASE::fi1] = - (curApprox[BASE::Ie] + C1 * curApprox[BASE::_UC1] + curApprox[BASE::Il]);
  free_part[BASE::fi2] = - ( - C1 * curApprox[BASE::UC1]
                  - curApprox[BASE::Il]
                  - It1 * (exp((curApprox[BASE::fi4] - curApprox[BASE::fi2]) / Mft) - 1.0)
                  - Cb1 * curApprox[BASE::_UCb1]
                  - (curApprox[BASE::fi4] - curApprox[BASE::fi2]) / Ru
                  + (curApprox[BASE::fi2] - curApprox[BASE::fi5]) / RB);
  free_part[BASE::fi3] = - ( - It2 * (exp((curApprox[BASE::fi5] - curApprox[BASE::fi3]) / Mft) - 1.0)
                  - Cb2 * curApprox[BASE::_UCb2]
                  - (curApprox[BASE::fi5] - curApprox[BASE::fi3]) / Ru
                  + C2 * curApprox[BASE::_UC2]
                  + curApprox[BASE::fi3] / R);
  free_part[BASE::fi4] = - ( curApprox[BASE::fi4] / RB
                  + It1 * (exp((curApprox[BASE::fi4] - curApprox[BASE::fi2]) / Mft) - 1.0)
                  + Cb1 * curApprox[BASE::_UCb1]
                  + (curApprox[BASE::fi4] - curApprox[BASE::fi2]) / Ru);
  free_part[BASE::fi5] = - (- (curApprox[BASE::fi2] - curApprox[BASE::fi5]) / RB
                  + It2 * (exp((curApprox[BASE::fi5] - curApprox[BASE::fi3]) / Mft) - 1.0)
                  + Cb2 * curApprox[BASE::_UCb2]
                  + (curApprox[BASE::fi5] - curApprox[BASE::fi3]) / Ru);
  free_part[BASE::Ie] = - (curApprox[BASE::fi1] - E(t));
}

void Scheme::show_global_matx(){
  for(int i = 0; i < size; i++){
    for(int j = 0; j < size; j++)
      cout << global_matx[i][j] << "\t\t";
    cout << endl << endl;
  }
}

void Scheme::show_free_part(){
  for(int i = 0; i < size; i++)
    cout << free_part[i] << endl;
}

void Scheme::complete_iter(double *x, double *approx, double *beginapprox){
  for(int i = 0; i < N; i++)
    approx[i] += x[i];
}

void complete_step(double *approx, double *beginapprox){
  for(int i = 0; i < N; i++){
    cout << "prev == " << beginapprox[i] << " " << " new  == " << approx[i] << endl;
  }
  for(int i = 0; i < N; i++){
    U_m_2[i] = U_m_1[i];
    U_m_1[i] = beginapprox[i];
    beginapprox[i] = approx[i];
  }
}

double* gausse(double **a, double *y, int n) {
	double *x, max;
  int k, index;
  const double eps = 1e-32;  // точность
  x = new double[n];
  k = 0;
  while (k < n) {
    max = fabs(a[k][k]);
    index = k;
  	for (int i = k + 1; i < n; i++) {
      if (fabs(a[i][k]) > max){
        max = fabs(a[i][k]);
        index = i;
      }
    }
    if (max < eps) {
      cout << "Решение получить невозможно из-за нулевого столбца ";
      cout << index << " матрицы A" << endl;
      return 0;
    }
    for (int j = 0; j < n; j++) {
      double temp = a[k][j];
      a[k][j] = a[index][j];
      a[index][j] = temp;
    }
    double temp = y[k];
    y[k] = y[index];
    y[index] = temp;
    for (int i = k; i < n; i++) {
  		double temp = a[i][k];
      if (fabs(temp) < eps)  continue; // для нулевого коэффициента пропустить
      for (int j = 0; j < n; j++)
        a[i][j] = a[i][j] / temp;
      	y[i] = y[i] / temp;
      	if (i == k)  continue; // уравнение не вычитать само из себя
      	for (int j = 0; j < n; j++)
        	a[i][j] = a[i][j] - a[k][j];
      	y[i] = y[i] - y[k];
    	}
    	k++;
    }
  for (k = n - 1; k >= 0; k--){
    x[k] = y[k];
    for (int i = 0; i < k; i++)
      y[i] = y[i] - a[i][k] * x[k];
  }
  return x;
}

int check_delta(double *x){
  int flag = 0;
  for(int i = 0; i < N; i++){
    if(fabs(x[i]) < delta)
      flag++;
  }
  return (flag == N) ? 1 : 0;
}

double search_min(double *x){
  double min = x[0];
  for(int i = 1; i < N; i++){
    if(min > x[i])
      min = x[i];
  }
  return min;
}

void take_solution1(){
  dt = dt / 2;
  for(int i = 0; i < N; i++)
    curApprox[i] = 0.0;
}

void take_solution2(){
  double cur_eps;
  int kek;
  for(int i = 0; i < N; i++){
    eps[BASE::fi3] = fabs( dt / (dt + dt_prev) * ((curApprox[BASE::fi3] - U_m_1[BASE::fi3]) - dt / dt_prev * (U_m_1[BASE::fi3] - U_m_2[BASE::fi3])));
  }
  cur_eps = search_min(eps);
  cur_eps = eps[BASE::fi3];
  cout << cur_eps << " CUR EPS" << endl;
  cin >> kek;
  if(cur_eps < EPS1){
    dt = dt * 2;
    complete_step(curApprox, prevStep);
    t += dt;
    dt_prev = dt;
  }
  if(cur_eps < EPS2){
    complete_step(curApprox, prevStep);
    t += dt;
    dt_prev = dt;
  }else{
    dt = dt / 2;
  }
}

int main(){
  double *x = new double[N];
  FILE *data;
  int end = 0;
  int iter_step = 0;
  int t_step = 0;
  int kek;
  bool is_continue = false;
  double time = 0.0;
  if((data = fopen("data.txt", "w"))==NULL){
    cout << "Can't open file" << endl;
  }

  //curApprox[BASE::_UCb1] = -6.0;
  //curApprox[BASE::_UCb2] = 6.0;

  //prevStep[BASE::_UCb1] = -6.0;
  //prevStep[BASE::_UCb2] = 6.0;

  while(t < TIME){
    while(!end){
      Scheme *scheme = new Scheme(N);
      scheme->init_global_matx();
      scheme->init_free_part();
      //scheme->show_global_matx();
      //scheme->show_free_part();
      if(iter_step > MAX_ITER_STEPS){
        cin >> kek;
        take_solution1();
        iter_step = 0;
        continue;
      }
      x = gausse(scheme->get_global_matx(), scheme->get_free_part(), N);
      cout << endl << "delta" << endl;
      for(int i = 0; i < N; i++)
        cout << x[i] << endl;
      for(int i = 0; i < N; i++){
        cout << "cur " << str[i] << " " << curApprox[i];
        cout << "\t" << "new " << curApprox[i] << " + " << x[i] << " == "<< curApprox[i] + x[i] << endl;
      }
      scheme->complete_iter(x, curApprox, prevStep);
      end = check_delta(x);
      iter_step++;
      cout << "STEP == " << iter_step << endl;
      cout << "TIME STEP == " << t_step << endl;
      cout << "TIME == " << t << endl;
      //cin >> kek;
    }
    if(t_step != 0)
      fprintf(data, "%e\n", curApprox[BASE::fi3] );
    t_step++;
    iter_step = 0;
    end = 0;
    for(int i = 0; i < N; i++){
      if(fabs(curApprox[BASE::fi3] - prevStep[BASE::fi3]) > INTEGR_DELTA/*7*/ && t_step - 1 != 1){
        cout << "solution 2" << endl;
        cin >> kek;
        take_solution2();
        is_continue = true;
        break;
      }
    }
    if(is_continue){
      is_continue = false;
      continue;
    }
    complete_step(curApprox, prevStep);
    dt_prev = dt;
    t += dt;
  }
  fclose(data);
  return 0;
}

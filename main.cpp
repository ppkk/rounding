#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <limits>
#include <typeinfo>

#define x2 ((x)*(x))
#define x3 ((x)*(x)*(x))
#define x4 ((x)*(x)*(x)*(x))

#define sqr(x) ((x) * (x))

using namespace std;

typedef float real_assemble;
typedef float real_solve;
typedef float real_error;


const real_assemble eps = 0;

real_error exact(real_error x)
{
   return x * (1-x);
}

real_error exact_der(real_error x)
{
   return 1 - 2 * x;
}

// used in assemble
real_assemble primitive_left(real_assemble a, real_assemble x, real_assemble h)
{
   return (-eps/4 * x4 + eps*(a+1)/3 * x3 + (2-eps*a)/2 * x2 - 2*a*x)/h;
}

real_assemble primitive_right(real_assemble a, real_assemble x, real_assemble h)
{
   return (eps/4 * x4 - eps*(a+h+1)/3 * x3 + (eps*(a+h)-2)/2 * x2 + 2*(a+h)*x)/h;
}

real_assemble rhs_val(int idx, real_assemble h)
{
   real_assemble a = (idx+1) * h;
   real_assemble ret = 0;
   ret += primitive_left(a-h, a, h) - primitive_left(a-h, a-h, h);
   ret += primitive_right(a, a+h, h) - primitive_right(a, a, h);
   return ret;
}

// gauss quadrature is used in error calculation
const int GAUSS_NP = 3;
real_error weights[GAUSS_NP] = {5./9., 8./9., 5./9.};
real_error gauss_points[GAUSS_NP] = {-sqrt(3./5.), 0, sqrt(3./5.)};


void plot_error(int N, real_error* sol, ofstream &f_error)
{
   const int intervals = N+1;
   const real_error h = 1./intervals;

   for (int interval_idx = 0; interval_idx < intervals; interval_idx++)
   {
      real_error a, b, a_val, b_val;
      a = h*interval_idx;
      b = h*(interval_idx+1);

      if(interval_idx == 0)
         a_val = 0.;
      else
         a_val = sol[interval_idx - 1];
      if(interval_idx == intervals - 1)
         b_val = 0.;
      else
         b_val = sol[interval_idx];

      for(real_error pt = -1; pt < 1; pt += 0.05)
      {
         real_error trans_pt = (b-a)/2 * pt + (a+b)/2;
         real_error trans_val = (b_val-a_val)/2 * pt + (a_val+b_val)/2;
         f_error << trans_pt << " " << exact(trans_pt) - trans_val << endl;
      }

   }
}

void add_error_contribution(real_error a, real_error b, real_error a_val, real_error b_val, real_error &l2, real_error &h1_semi, ofstream &f_error)
{
   real_error derivative = (b_val - a_val) / (b-a);

   for (int quad_idx = 0; quad_idx < GAUSS_NP; quad_idx++)
   {
      real_error trans_pt = (b-a)/2 * gauss_points[quad_idx] + (a+b)/2;
      real_error trans_val = (b_val-a_val)/2 * gauss_points[quad_idx] + (a_val+b_val)/2;

      l2 += (b-a)/2 * weights[quad_idx] * sqr(exact(trans_pt) - trans_val);
      h1_semi += (b-a)/2 * weights[quad_idx] * sqr(exact_der(trans_pt) - derivative);

      //f_error << trans_pt << " " << exact(trans_pt) - trans_val << endl;
   }
}

void calc_errors(int N, real_error* sol, real_error &error_l2, real_error &error_h1_semi)
{
   string str_sol = "sol/solution-" + to_string(N) + ".txt";
   ofstream f_solution(str_sol);
   string str_err = "sol/error-" + to_string(N) + ".txt";
   ofstream f_error(str_err);
   const int intervals = N+1;
   const real_error h = 1./intervals;

   error_l2 = 0.0;
   error_h1_semi = 0.0;

   real_error left_pt, right_pt, left_val, right_val;
   for(int interval_idx = 0; interval_idx < intervals; interval_idx++)
   {
      left_pt = h*interval_idx;
      right_pt = h*(interval_idx+1);

      if(interval_idx == 0)
         left_val = 0.;
      else
         left_val = sol[interval_idx - 1];
      if(interval_idx == intervals - 1)
         right_val = 0.;
      else
         right_val = sol[interval_idx];

      add_error_contribution(left_pt, right_pt, left_val, right_val, error_l2, error_h1_semi, f_error);
   }

   error_l2 = sqrt(error_l2);
   error_h1_semi = sqrt(error_h1_semi);


//   plot_error(N, sol, f_error);
//   f_solution <<std::setprecision (std::numeric_limits<real_error>::digits10) << 0.0 << " " << 0.0 << endl;
//   for(int i = 0; i < N; i++)
//      f_solution << h*(i+1) << " " << sol[i] << endl;
//   f_solution << 1.0 << " " << 0.0 << endl;
   f_solution.close();
   f_error.close();
}

void assemble(int N, real_assemble *diag, real_assemble *upper, real_assemble *lower, real_assemble* rhs)
{

   const int intervals = N+1;
   const real_assemble h = 1./intervals;

   real_assemble u_v_same = h/3.;
   real_assemble u_v_diff = h/6.;
   real_assemble du_dv_same = 1./h;
   real_assemble du_dv_diff = -1./h;

   //cout <<std::setprecision (std::numeric_limits<real_assemble>::digits10) << "h " << h << ", uv " << u_v_same << ", duv "<< du_dv_same << ", eps*uv + duv " << eps*u_v_same + du_dv_same << endl;

   for(int i = 0; i < N; i++)
      diag[i] = 2 * (eps * u_v_same + du_dv_same);

   for(int i = 0; i < N - 1; i++)
   {
      upper[i] = eps * u_v_diff + du_dv_diff;
      lower[i] = eps * u_v_diff + du_dv_diff;
   }

   for(int i = 0; i < N; i++)
   {
      rhs[i] = rhs_val(i, h);
      //cout << "rhs: " << rhs[i] << endl;
   }
}


void solve(int N, real_solve *diag, real_solve *upper, real_solve *lower, real_solve* rhs, real_solve* sol)
{
   for(int i = 0; i < N; i++)
   {
      sol[i] = 0;
   }

   // forward elimination:
   for(int i = 0; i < N-1; i++)
   {
      real_solve multiple = -lower[i]/diag[i];
      diag[i+1] += multiple * upper[i];
      lower[i] += multiple * diag[i];
      rhs[i+1] += multiple * rhs[i];
   }

//   for(int i = 0; i < N-1; i++)
//      assert(fabs(lower[i]) < 1e-10);

   // backward elimination:
   sol[N-1] = rhs[N-1] / diag[N-1];
   for(int i = N-2; i >= 0; i--)
   {
      sol[i] = (rhs[i] - upper[i] * sol[i+1]) / diag[i];
   }

}

void run(int N, real_error &error_l2, real_error &error_h1_semi)
{
   real_assemble* diag_assemble = new real_assemble[N];
   real_assemble* upper_assemble = new real_assemble[N - 1];
   real_assemble* lower_assemble = new real_assemble[N - 1];
   real_assemble* rhs_assemble = new real_assemble[N];

   assemble(N, diag_assemble, upper_assemble, lower_assemble, rhs_assemble);

   real_solve* diag_solve = new real_solve[N];
   real_solve* rhs_solve = new real_solve[N];
   for(int i = 0; i < N; i++)
   {
      diag_solve[i] = (real_solve)diag_assemble[i];
      rhs_solve[i] = (real_solve)rhs_assemble[i];
   }
   delete[] diag_assemble;
   delete[] rhs_assemble;

   real_solve* upper_solve = new real_solve[N - 1];
   real_solve* lower_solve = new real_solve[N - 1];
   for(int i = 0; i < N - 1; i++)
   {
      upper_solve[i] = (real_solve)upper_assemble[i];
      lower_solve[i] = (real_solve)lower_assemble[i];
   }
   delete[] upper_assemble;
   delete[] lower_assemble;

   real_solve* sol_solve = new real_solve[N];

   solve(N, diag_solve, upper_solve, lower_solve, rhs_solve, sol_solve);

   delete[] diag_solve;
   delete[] upper_solve;
   delete[] lower_solve;
   delete[] rhs_solve;

   real_error* sol_error = new real_error[N];
   for(int i = 0; i < N; i++)
   {
      sol_error[i] = (real_error)sol_solve[i];
   }
   delete[] sol_solve;

   calc_errors(N, sol_error, error_l2, error_h1_semi);

   delete[] sol_error;
}

int main()
{
//   cout << sizeof(real_assemble) << ", " << sizeof(real_solve) << ", "<< sizeof(real_error) << endl;

   for(double N = 2; N < 1e7; N = round(N*1.3))
   {
      real_error error_l2, error_h1_semi;
      run(N, error_l2, error_h1_semi);
//      cout << "N: " << N << ", error " << error << endl;
      cout << N << " " << error_l2 << " " << error_h1_semi << endl;
   }

   return 0;
}

int main2()
{
   //cout <<  std::numeric_limits<double>::digits10 << endl;
   real_error error_l2, error_h1_semi;
   run(300, error_l2, error_h1_semi);
   return 0;
}

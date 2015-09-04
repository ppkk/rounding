#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <limits>
#include <typeinfo>

#include "quadrature.h"
#include "mesh.h"

/// equation -u'' + eps u = f(x)
/// exact solution u(x) = x*(1-x) = x - x^2
/// rhs 2 + eps (x - x^2)

#define sqr(x) ((x) * (x))

using namespace std;

typedef float real_assemble;
typedef double real_solve;
typedef float real_error;


const real_assemble eps = 1;

real_error exact(real_error x)
{
   return x * (1-x);
}

real_error exact_der(real_error x)
{
   return 1 - 2 * x;
}

real_assemble rhs_val(real_assemble x)
{
   return (real_assemble(eps) * (x*(1-x))) + 2;
}


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

void add_error_contribution(real_error a, real_error b, real_error a_val, real_error b_val, real_error &l2, real_error &h1_semi, ofstream &f_error, Quadrature<real_error> &quad)
{
   real_error derivative = (b_val - a_val) / (b-a);

   for (int quad_idx = 0; quad_idx < quad.n_points; quad_idx++)
   {
      real_error trans_pt = (b-a)/2 * quad.points[quad_idx] + (a+b)/2;
      real_error trans_val = (b_val-a_val)/2 * quad.points[quad_idx] + (a_val+b_val)/2;

      l2 += (b-a)/2 * quad.weights[quad_idx] * sqr(exact(trans_pt) - trans_val);
      h1_semi += (b-a)/2 * quad.weights[quad_idx] * sqr(exact_der(trans_pt) - derivative);

//      f_error << trans_pt << " " << exact(trans_pt) - trans_val << endl;
   }
}

void calc_errors(int N, real_error* sol, real_error &error_l2, real_error &error_h1_semi)
{
   Quadrature<real_error> quad(3);

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

      add_error_contribution(left_pt, right_pt, left_val, right_val, error_l2, error_h1_semi, f_error, quad);
   }

   error_l2 = sqrt(error_l2);
   error_h1_semi = sqrt(error_h1_semi);


   plot_error(N, sol, f_error);
   f_solution <<std::setprecision (std::numeric_limits<real_error>::digits10) << 0.0 << " " << 0.0 << endl;
   for(int i = 0; i < N; i++)
      f_solution << h*(i+1) << " " << sol[i] << endl;
   f_solution << 1.0 << " " << 0.0 << endl;
   f_solution.close();
   f_error.close();
}


void assemble(Mesh<real_assemble> &mesh, real_assemble *diag, real_assemble *upper, real_assemble *lower, real_assemble* rhs)
{
//   cout <<std::setprecision (std::numeric_limits<real_assemble>::digits10) << "h " << mesh.h << ", uv " << u_v_same << ", duv "<< du_dv_same << ", eps*uv + duv " << eps*u_v_same + du_dv_same << endl;

   for(int i = 0; i < mesh.N; i++)
   {
      diag[i] = 0.0;
      rhs[i] = 0.0;
   }

   for(int i = 0; i < mesh.N - 1; i++)
   {
      upper[i] = 0.0;
      lower[i] = 0.0;
   }

   Quadrature<real_assemble> quad(3);

   for(int interval_idx = 0; interval_idx < mesh.num_intervals; interval_idx++)
   {
      real_assemble a, b;
      mesh.interval_endpoints(interval_idx, a, b);
      real_assemble h = b - a;

      real_assemble u_v_same = h/3.;
      real_assemble u_v_diff = h/6.;
      real_assemble du_dv_same = 1./h;
      real_assemble du_dv_diff = -1./h;

      real_assemble diag_contrib = eps * u_v_same + du_dv_same;
      real_assemble offdiag_contrib = eps * u_v_diff + du_dv_diff;

      if(interval_idx > 0)
      {
         diag[interval_idx - 1] += diag_contrib;
      }
      if(interval_idx < mesh.num_intervals - 1)
      {
         diag[interval_idx] += diag_contrib;
      }
      if((interval_idx > 0) && (interval_idx < mesh.num_intervals - 1))
      {
         upper[interval_idx - 1] += offdiag_contrib;
         lower[interval_idx - 1] += offdiag_contrib;
      }

      for (int quad_idx = 0; quad_idx < quad.n_points; quad_idx++)
      {
         real_error trans_pt = (b-a)/2 * quad.points[quad_idx] + (a+b)/2;
         real_error trans_val_left_basis = (0.-1.)/2 * quad.points[quad_idx] + (1.+0.)/2;
         real_error trans_val_right_basis = (1.-0.)/2 * quad.points[quad_idx] + (0.+1.)/2;
         if(interval_idx > 0)
         {
            rhs[interval_idx - 1] += (b-a)/2 * quad.weights[quad_idx] * trans_val_left_basis * rhs_val(trans_pt);
         }
         if(interval_idx < mesh.num_intervals - 1)
         {
            rhs[interval_idx] += (b-a)/2 * quad.weights[quad_idx] * trans_val_right_basis * rhs_val(trans_pt);
         }
      }
   }

//   for(int i = 0; i < N; i++)
//   {
//      cout << "rhs " << rhs[i] << endl;
//   }

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

   Mesh<real_assemble> mesh_assemble;
   mesh_assemble.create_uniform(N+1);

   assemble(mesh_assemble, diag_assemble, upper_assemble, lower_assemble, rhs_assemble);

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

   Mesh<real_solve> mesh_solve;
   convert_mesh(mesh_assemble, mesh_solve);
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
   Mesh<real_assemble> mesh_assemble;
   //cout <<  std::numeric_limits<double>::digits10 << endl;
   real_error error_l2, error_h1_semi;
   //run(30, error_l2, error_h1_semi);


   return 0;
}

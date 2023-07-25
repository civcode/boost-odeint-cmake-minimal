#include <cmath>
#include <iostream>
#include <iomanip>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#include "tic_toc_timer.h"

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

//typedef boost::array<double , 3> state_type;
typedef boost::array<double , 5> state_type;
typedef runge_kutta_dopri5< double > stepper_type;
typedef runge_kutta4<state_type> rk4;
typedef runge_kutta_cash_karp54< state_type > rkck54;
typedef controlled_runge_kutta< rkck54 > ctrl_rkck54;

using std::cout;
using std::endl;


void lorenz(const state_type &x , state_type &dxdt , double t)
{
  dxdt[0] = sigma * (x[1] - x[0]);
  dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
  dxdt[2] = -b * x[2] + x[0] * x[1];
}

void model(const state_type &x, state_type &dxdt, double t) 
{
  double L1 = 0.5;
  double L2 = 0.5;
  
  double v = 1.00;
  double wp = 0.0;

  double theta = x[2];
  double gamma = x[4];

  double c2 = (L1 * std::cos(gamma) + L2);

  dxdt[0] = v * std::cos(theta);
  dxdt[1] = v * std::sin(theta);
  dxdt[2] = (v * std::sin(gamma) / c2) + (wp * L2/c2);
  dxdt[3] = (v * std::sin(gamma) / c2) - (wp * L1*std::cos(gamma)/c2);
  dxdt[4] = wp;

}

void rhs( const double x , double &dxdt , const double t )
{
    dxdt = 3.0/(2.0*t*t) + x/(2.0*t);
}

void write_lorenz(const state_type &x , const double t)
{
  //cout << std::setw(3);
  //cout << cout.precision(2); // 
  cout << std::setprecision(2);
  //cout << cout.width(2);
  //cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << " " << x[4] << endl;
  cout << t << '\t';
  for (auto elem : x) {
    //cout << t << '\t' << elem << '\t';
    cout << elem << '\t';
  }
  cout << endl;
}

int main(int argc, char **argv)
{
  TicTocTimer timer;

  //state_type x = {10.0 , 1.0 , 1.0}; // initial conditions
  state_type x = {0.0 , 0.0 , 0.0, 0.0, 0.76}; // initial conditions
  state_type y = {0.0 , 0.0 , 0.0, 0.0, 0.76}; // initial conditions
  size_t step_count;
  //double x = 0;
  //integrate(lorenz, x , 0.0 , 25.0 , 0.1 , write_lorenz);
  //integrate(lorenz, x, 0.0, 25.0, 0.1);
  cout << "t\t" << "x\t" << "y\t" << "th1\t" << "th2\t" << "gm\t" << endl;
  timer.tic();
  //size_t step_count = integrate(model, x, 0.0, 5.0, 0.01, write_lorenz);
  //size_t step_count = integrate_adaptive(make_controlled(1E-12, 1E-12, stepper_type()), rhs, x, 0.0, 5.0, 0.01);
  //size_t step_count = integrate_adaptive(rk4(), model, x, 0.0, 5.0, 0.1, write_lorenz);
  //size_t step_count = integrate_adaptive(rk4(), model, x, 0.0, 5.0, 0.01, write_lorenz);
  //size_t step_count = integrate_adaptive(rk4(), model, x, 0.0, 5.0, 0.01, write_lorenz);
  //size_t step_count = integrate_adaptive(ctrl_rkck54(), model, x, 0.0, 5.0, 0.01, write_lorenz);
  step_count = integrate_adaptive(rk4(), model, x, 0.0, 1.0, 0.001, write_lorenz);
  cout << "time expired: " << timer.toc().ms().to_string() << endl;
  cout << "step count:   " << step_count << endl;
  timer.tic();
  step_count = integrate_adaptive(ctrl_rkck54(), model, y, 0.0, 1.0, 0.001, write_lorenz);
  //cout << "time expired: " << timer.toc().ms().value() << endl;
  cout << "time expired: " << timer.toc().ms().to_string() << endl;
  cout << "step count:   " << step_count << endl;
}

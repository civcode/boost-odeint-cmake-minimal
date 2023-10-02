
#include <cmath>
#include <iostream>
#include <iomanip>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#include "tic_toc_timer.h"

using namespace std;
using namespace boost::numeric::odeint;
using namespace std::placeholders;

typedef std::vector< double > state_type;

/* The rhs of x' = f(x) defined as a class */
class harm_osc {

    double m_gam;

public:
    static double gam_;
    harm_osc( double gam ) : m_gam(gam) { }

    void operator() ( const state_type &x , state_type &dxdt , const double /* t */ )
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - m_gam*x[1];
    }
    static void func_static ( const state_type &x , state_type &dxdt , const double /* t */ )
    {
        double m_gam = 1.0;
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - m_gam*x[1];
    }

    void func_non_static ( const state_type &x , state_type &dxdt , const double /* t */ )
    {
        m_gam = 1.0;
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - m_gam*x[1];
    }

    // static void wrapper(const state_type &x , state_type &dxdt , const double t) {

    //   this->operator(x, dxdt, t);
    // }

    void run_static() {
      state_type x = {0.10 , 0.0}; // initial conditions
      state_type y = {0.10 , 0.0}; // initial conditions
      //size_t step_count = integrate(this->operator(), x, 0.0, 5.0, 0.01, write_lorenz);
      size_t step_count = integrate(this->func_static, x, 0.0, 5.0, 0.01, write_lorenz);
      //size_t step_count = integrate(std::bind(&harm_osc::func_non_static, this, _1, _2, _3), x, 0.0, 5.0, 0.01, write_lorenz);

    }

    void run_non_static() {
      state_type x = {0.10 , 0.0}; // initial conditions
      state_type y = {0.10 , 0.0}; // initial conditions
      //size_t step_count = integrate(this->operator(), x, 0.0, 5.0, 0.01, write_lorenz);
      //size_t step_count = integrate(this->func, x, 0.0, 5.0, 0.01, write_lorenz);
      size_t step_count = integrate(std::bind(&harm_osc::func_non_static, this, _1, _2, _3), x, 0.0, 5.0, 0.01, write_lorenz);

    }

    static void write_lorenz(const state_type &x , const double t)
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
};


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


double harm_osc::gam_ = 0.0;

int main(int argc, char **argv)
{
  TicTocTimer timer;

  harm_osc osc(1.0);
  harm_osc::gam_ = 1.0;

  //state_type x = {10.0 , 1.0 , 1.0}; // initial conditions
  //state_type x = {0.0 , 0.0 , 0.0, 0.0, 0.76}; // initial conditions
  //state_type y = {0.0 , 0.0 , 0.0, 0.0, 0.76}; // initial conditions
  state_type x = {0.0 , 0.0}; // initial conditions
  state_type y = {0.0 , 0.0}; // initial conditions
  //size_t step_count;
  //double x = 0;
  //integrate(lorenz, x , 0.0 , 25.0 , 0.1 , write_lorenz);
  //integrate(lorenz, x, 0.0, 25.0, 0.1);
  cout << "t\t" << "x\t" << "y\t" << "th1\t" << "th2\t" << "gm\t" << endl;
  timer.tic();
  //size_t step_count = integrate(model, x, 0.0, 5.0, 0.01, write_lorenz);
  size_t step_count = integrate(osc, x, 0.0, 5.0, 0.01, write_lorenz);
  //size_t step_count = integrate_adaptive(make_controlled(1E-12, 1E-12, stepper_type()), rhs, x, 0.0, 5.0, 0.01);
  //size_t step_count = integrate_adaptive(rk4(), model, x, 0.0, 5.0, 0.1, write_lorenz);
  //size_t step_count = integrate_adaptive(rk4(), model, x, 0.0, 5.0, 0.01, write_lorenz);
  //size_t step_count = integrate_adaptive(rk4(), model, x, 0.0, 5.0, 0.01, write_lorenz);
  //size_t step_count = integrate_adaptive(ctrl_rkck54(), model, x, 0.0, 5.0, 0.01, write_lorenz);
 // step_count = integrate_adaptive(rk4(), model, x, 0.0, 1.0, 0.001, write_lorenz);
  cout << "time expired: " << timer.toc().ms().to_string() << endl;
  cout << "step count:   " << step_count << endl;
  timer.tic();
  //step_count = integrate_adaptive(ctrl_rkck54(), model, y, 0.0, 1.0, 0.001, write_lorenz);
  //cout << "time expired: " << timer.toc().ms().value() << endl;
  cout << "time expired: " << timer.toc().ms().to_string() << endl;
  cout << "step count:   " << step_count << endl;

  cout << "static" << endl;
  timer.tic();
  osc.run_static();
  cout << "time expired: " << timer.toc().ms().to_string() << endl;

  cout << "non-static" << endl;
  timer.tic();
  osc.run_non_static();
  cout << "time expired: " << timer.toc().ms().to_string() << endl;
}

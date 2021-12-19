#include <chrono> // for high_resolution_clock
#include <ctime>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <iostream>
#include <opstop/ipopt_problem.hpp>

using namespace opstop;
int main() {
  std::size_t number_of_wp = 10;
  std::size_t codom_dim = 7;
  std::size_t nglp = 10;

  double acc_max = 10.0;

  double exec_time = (double)number_of_wp - 1.0;
  double ti = 0.5 * exec_time;

  Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));

  gsplines::GSpline trj = gsplines::optimization::optimal_sobolev_norm(
      wp, gsplines::basis::BasisLegendre(6), {{1.0, 3}}, exec_time);

  gsplines::functions::FunctionExpression diffeo =
      minimum_time_bouded_acceleration(trj, ti, 10);

  gsplines::functions::FunctionExpression diffeo_diff_1 = diffeo.derivate();
  gsplines::functions::FunctionExpression diffeo_diff_2 =
      diffeo_diff_1.derivate();

  double T = diffeo.get_domain().second;

  Eigen::VectorXd time_span = Eigen::VectorXd::LinSpaced(10, ti, T);

  Eigen::MatrixXd diffeo_hist = diffeo(time_span);

  Eigen::MatrixXd diffeo_diff_1_hist = diffeo_diff_1(time_span);
  Eigen::MatrixXd diffeo_diff_2_hist = diffeo_diff_2(time_span);

  Eigen::VectorXd vec(1);
  vec(0) = ti;

  double ti_test = diffeo(vec)(0, 0);
  double ve = diffeo_diff_1(vec)(0, 0);
  double ac = diffeo_diff_2(vec)(0, 0);

  diffeo.print();
  assert(std::abs(ti_test - ti) < 1.0e-9);
  assert(std::abs(ve - 1.0) < 1.0e-9);
  assert(std::abs(ac) < 1.0e-9);

  gsplines::functions::FunctionExpression stop_trj = trj.compose(diffeo);
  gsplines::functions::FunctionExpression stop_trj_diff_1 = stop_trj.derivate();
  gsplines::functions::FunctionExpression stop_trj_diff_2 =
      stop_trj.derivate(2);

  vec(0) = T;
  ti_test = diffeo(vec)(0, 0);
  ve = diffeo_diff_1(vec)(0, 0);
  ac = diffeo_diff_2(vec)(0, 0);

  assert(std::abs(ve) < 1.0e-9);
  assert(std::abs(ac) < 1.0e-9);

  double vel_value = stop_trj_diff_1(vec).row(0).norm();
  double acc_value = stop_trj_diff_2(vec).row(0).norm();

  assert(vel_value < 1.0e-9);
  assert(acc_value < 1.0e-9);

  double mean_time = 0.0;
  std::size_t number_of_tests = 1;
  double time_elapsed_ms = 0;
  for (std::size_t _ = 0; _ < number_of_tests; _++) {
    std::clock_t t_cpu_1 = std::clock();
    auto t1 = std::chrono::high_resolution_clock::now();
    gsplines::functions::FunctionExpression _diffeo =
        minimum_time_bouded_acceleration(trj, ti, acc_max);
    gsplines::functions::FunctionExpression _diffeo_diff_1 = _diffeo.derivate();
    gsplines::functions::FunctionExpression _diffeo_diff_2 =
        _diffeo_diff_1.derivate();

    double _T = diffeo.get_domain().second;

    Eigen::VectorXd _time_span =
        Eigen::VectorXd::LinSpaced((T - ti) / 0.01, ti, T);

    Eigen::MatrixXd _diffeo_hist = diffeo(time_span);

    Eigen::MatrixXd _diffeo_diff_1_hist = diffeo_diff_1(time_span);
    Eigen::MatrixXd _diffeo_diff_2_hist = diffeo_diff_2(time_span);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::clock_t t_cpu_2 = std::clock();
    time_elapsed_ms = 1000.0 * (t_cpu_2 - t_cpu_1) / CLOCKS_PER_SEC;
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    mean_time += ms_double.count();
  }
  mean_time /= number_of_tests;
  assert(mean_time < 100.0);
  printf("Mean optimization time is %+14.7e ms\n", mean_time);
  printf("Mean optimization time is %+14.7e ms CPU\n", time_elapsed_ms);
  return 0;
}

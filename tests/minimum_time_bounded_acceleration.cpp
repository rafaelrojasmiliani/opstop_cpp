#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/ipopt_solver.hpp>
#include <iostream>
#include <opstop/ipopt_problem.hpp>

using namespace opstop;
int main() {
  std::size_t number_of_wp = 3;
  std::size_t codom_dim = 7;
  std::size_t nglp = 10;

  double exec_time = (double)number_of_wp - 1.0;
  double ti = 0.5 * exec_time;

  Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));

  gsplines::PiecewiseFunction trj = gsplines_opt::optimal_sobolev_norm(
      wp, gsplines::basis::BasisLegendre(6), {{4, 1}}, exec_time);

  gsplines::functions::FunctionExpression diffeo =
      minimum_time_bouded_acceleration(trj, ti, 10);

  gsplines::functions::FunctionExpression diffeo_diff_1 = diffeo.derivate();
  gsplines::functions::FunctionExpression diffeo_diff_2 =
      diffeo_diff_1.derivate();

  double T = diffeo.get_domain().second;

  Eigen::VectorXd time_span = Eigen::VectorXd::LinSpaced(10, ti, T);

  Eigen::MatrixXd diffeo_hist = diffeo(time_span);

  std::cout << "after the call \n----\n" << diffeo_hist << "\n----------\n";
  Eigen::MatrixXd diffeo_diff_1_hist = diffeo_diff_1(time_span);
  Eigen::MatrixXd diffeo_diff_2_hist = diffeo_diff_2(time_span);

  Eigen::VectorXd vec(1);
  vec(0) = ti;

  double ve = diffeo_diff_1(vec)(0, 0);
  double ac = diffeo_diff_2(vec)(0, 0);

  assert(std::abs(diffeo(vec)(0, 0) - ti) < 1.0e-9);
  assert(std::abs(diffeo_diff_1(vec)(0, 0) - 1.0) < 1.0e-9);
  assert(std::abs(diffeo_diff_2(vec)(0, 0)) < 1.0e-9);

  gsplines::functions::FunctionExpression stop_trj = trj.compose(diffeo);
  gsplines::functions::FunctionExpression stop_trj_diff_1 = trj.derivate();
  gsplines::functions::FunctionExpression stop_trj_diff_2 = trj.derivate(2);

  vec(0) = T;

  assert(std::abs(diffeo_diff_1(vec)(0, 0)) < 1.0e-9);
  assert(std::abs(diffeo_diff_2(vec)(0, 0)) < 1.0e-9);

  Eigen::MatrixXd vel_value = stop_trj_diff_1(vec);
  Eigen::MatrixXd acc_value = stop_trj_diff_2(vec);

  assert(vel_value.row(0).norm() < 1.0e-9);
  assert(acc_value.row(0).norm() < 1.0e-9);

  return 0;
}

#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/ipopt_solver.hpp>
#include <opstop/ipopt_problem.hpp>

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

  double T = diffeo.get_domain().second;

  gsplines::functions::FunctionExpression stop_trj = trj.compose(diffeo);
  gsplines::functions::FunctionExpression stop_trj_diff_1 = trj.derivate();
  gsplines::functions::FunctionExpression stop_trj_diff_2 = trj.derivate(2);

  Eigen::VectorXd vec(1);
  vec(0) = T;

  Eigen::MatrixXd vel_value = stop_trj_diff_1(vec);
  Eigen::MatrixXd acc_value = stop_trj_diff_2(vec);

  assert(vel_value.row(0).norm() < 1.0e-9);
  assert(acc_value.row(0).norm() < 1.0e-9);

  return 0;
}

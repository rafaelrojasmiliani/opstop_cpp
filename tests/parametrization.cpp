#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <gsplines++/ipopt_solver.hpp>
#include <opstop/parametrization.hpp>

std::size_t number_of_wp = 3;
std::size_t codom_dim = 3;
double exec_time = number_of_wp - 1;
Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));
gsplines::PiecewiseFunction trj = gsplines_opt::optimal_sobolev_norm(
    wp, gsplines::basis::BasisLegendre(6), {{4, 1}}, exec_time);

void test_pol_evaluation() {
  double Ts = 1;
  double ti = 0;
  double sf = 1;

  Eigen::VectorXd pol_coeff(6);

  pol_coeff << ti, Ts, 0.0, -6.0 * Ts + 10.0 * sf - 10.0 * ti,
      8.0 * Ts - 15.0 * sf + 15 * ti, -3.0 * Ts + 6.0 * sf - 6.0 * ti;

  gsplines::functions::CanonicPolynomial pol({0, 1}, pol_coeff);

  ParametrizedCurveHelper helper(trj, 9, ti);

  Eigen::VectorXd s_val = pol(helper.glp_);
  Eigen::VectorXd s_diff_1_wrt_tau = pol.derivate()(helper.glp_);

  helper.set_diffeo(Ts, sf);

  helper.compute_s_and_its_derivatives_wrt_tau();

  assert((helper.s_val_buff_ - pol(helper.glp_)).norm() < 1.0e-9);
  assert((helper.s_diff_1_wrt_tau_buff_ - pol.derivate()(helper.glp_)).norm() <
         1.0e-9);
  assert((helper.s_diff_2_wrt_tau_buff_ - pol.derivate(2)(helper.glp_)).norm() <
         1.0e-9);
  assert((helper.s_diff_3_wrt_tau_buff_ - pol.derivate(3)(helper.glp_)).norm() <
         1.0e-9);
}

int main() {

  test_pol_evaluation();

  return 0;
}

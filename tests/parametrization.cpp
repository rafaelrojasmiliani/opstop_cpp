#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <gsplines++/ipopt_solver.hpp>
#include <opstop/parametrization.hpp>

std::size_t number_of_wp = 3;
std::size_t codom_dim = 3;

double exec_time = (double)number_of_wp - 1.0;
double T = exec_time * 0.9;
double ti = 0.5 * exec_time;
double Ts = T - ti;
double sf = ti * (Ts / ti + 1);

Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));

gsplines::PiecewiseFunction trj = gsplines_opt::optimal_sobolev_norm(
    wp, gsplines::basis::BasisLegendre(6), {{4, 1}}, exec_time);

Eigen::VectorXd pol_coeff(6);

ParametrizedCurveHelper helper(trj, 9, ti);

void compare_assert(Eigen::MatrixXd &_m_nom, Eigen::MatrixXd &_m_test) {

  std::cout << _m_nom << "\n----\n";
  std::cout << _m_test << "\n----\n";
  std::cout << (_m_nom - _m_test).rowwise().norm() << "\n----\n";
  if (_m_nom.array().abs().maxCoeff() < 1.0e-9) {
    assert((_m_nom - _m_test).rowwise().norm().maxCoeff() < 1.0e-9);
  } else {
    double err = (_m_nom - _m_test).rowwise().norm().maxCoeff() /
                 _m_nom.rowwise().norm().maxCoeff();

    assert(err < 1.0e-9);
  }
}

void init() {
  pol_coeff << ti, Ts, 0.0, -6.0 * Ts + 10.0 * sf - 10.0 * ti,
      8.0 * Ts - 15.0 * sf + 15.0 * ti, -3.0 * Ts + 6.0 * sf - 6.0 * ti;
}

void test_pol_evaluation() {
  gsplines::functions::CanonicPolynomial pol({0, 1}, pol_coeff);

  Eigen::VectorXd s_val = pol(helper.glp_);
  Eigen::VectorXd s_diff_1_wrt_tau = pol.derivate()(helper.glp_);

  helper.set_diffeo(Ts, sf);
  helper.compute_s_and_its_derivatives_wrt_tau();

  assert((helper.s_val_buff_ - s_val).norm() < 1.0e-9);
  assert((helper.s_diff_1_wrt_tau_buff_ - s_diff_1_wrt_tau).norm() < 1.0e-9);
  assert((helper.s_diff_2_wrt_tau_buff_ - pol.derivate(2)(helper.glp_)).norm() <
         1.0e-9);
  assert((helper.s_diff_3_wrt_tau_buff_ - pol.derivate(3)(helper.glp_)).norm() <
         1.0e-9);
}

void test_compoisition_eval() {

  gsplines::functions::CanonicPolynomial pol({ti, exec_time}, pol_coeff);
  gsplines::functions::FunctionExpression tau_par =
      (1.0 / Ts) *
      (gsplines::functions::Identity({ti, exec_time}) +
       gsplines::functions::ConstFunction({ti, exec_time}, 1, -ti));
  gsplines::functions::FunctionExpression diffeo =
      gsplines::functions::Identity({0, ti}).concat(pol.compose(tau_par));
  gsplines::functions::FunctionExpression tau_par_inv =
      Ts * gsplines::functions::Identity({0, 1}) +
      gsplines::functions::ConstFunction({0, 1}, 1, ti);

  gsplines::functions::FunctionExpression exp = trj.compose(diffeo);
  gsplines::functions::FunctionExpression exp_2 = exp.compose(tau_par_inv);
  gsplines::functions::FunctionExpression exp_d1 =
      exp.derivate().compose(tau_par_inv);
  gsplines::functions::FunctionExpression exp_d2 =
      exp.derivate(2).compose(tau_par_inv);
  gsplines::functions::FunctionExpression exp_d3 =
      exp.derivate(3).compose(tau_par_inv);

  helper.set_diffeo(Ts, sf);
  helper.compute_s_and_its_derivatives_wrt_tau();
  helper.compute_q_and_its_derivatives_wrt_s();
  helper.compute_q_and_its_derivatives_wrt_t();

  Eigen::MatrixXd q_val = exp_2(helper.glp_);
  Eigen::MatrixXd q_d1 = exp_d1(helper.glp_);
  Eigen::MatrixXd q_d2 = exp_d2(helper.glp_);
  Eigen::MatrixXd q_d3 = exp_d3(helper.glp_);
  assert((helper.q_val_buff_ - q_val).norm() < 1.0e-9);

  assert((helper.q_diff_1_wrt_t_buff_ - q_d1).norm() < 1.0e-9);
  assert((helper.q_diff_2_wrt_t_buff_ - q_d2).norm() < 1.0e-9);

  compare_assert(helper.q_diff_3_wrt_t_buff_, q_d3);
}

int main() {

  init();

  test_pol_evaluation();

  test_compoisition_eval();

  return 0;
}

#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines++/ipopt_solver.hpp>
#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>
#include <opstop/torque_constraint.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include <opstop/differ.hpp>

#include <fenv.h>

using namespace opstop;
std::size_t number_of_wp = 3;
std::size_t codom_dim = 7;
std::size_t nglp = 10;

double exec_time = (double)number_of_wp - 1.0;
double T = exec_time * 0.9;
double ti = 0.5 * exec_time;
double Ts = 0.8 * T - ti;
double sf = ti * (Ts / ti + 1);

double eta_0 = 1.0 / 15.0 * (8.0 * exec_time / ti - 7.0);
double xi_0 = (exec_time / ti - 1.0);

double sf_radius = (exec_time / ti - eta_0) * ti / 2;
double sf_center = ((exec_time / ti - eta_0) / 2.0 + eta_0) * ti;
double Ts_center = xi_0 / 2.0 * ti;
double Ts_radius = xi_0 / 2.0 * ti;

Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));

gsplines::PiecewiseFunction trj = gsplines_opt::optimal_sobolev_norm(
    wp, gsplines::basis::BasisLegendre(6), {{4, 1}}, exec_time);

Eigen::VectorXd pol_coeff(6);

int main() {
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  pinocchio::Model model;

  pinocchio::urdf::buildModel("urdf/panda_arm.urdf", model);

  std::vector<double> tb(7, 10.0);
  TorqueConstraint cnstrt(trj, nglp, ti, tb, model);

  Eigen::VectorXd glp(10);

  std::tie(glp, std::ignore) =
      gsplines::collocation::legendre_gauss_lobatto_points_and_weights(10);

  Eigen::VectorXd sf_points =
      (0.5 * (eta_0 * ti + exec_time) + sf_radius * glp.array()).matrix();

  Eigen::VectorXd Ts_points =
      (0.5 * (xi_0 * ti) + Ts_radius * glp.array()).matrix();

  Eigen::Vector2d x;
  x(0) = Ts_center;
  x(1) = sf_center;
  Eigen::VectorXd tau = cnstrt.__GetValues(x);

  Eigen::MatrixXd jac_test = cnstrt.__GetJacobian(x);

  Eigen::MatrixXd jac_nom(jac_test);
  jac_nom.setConstant(0.0);

  std::size_t precision_order = 3;
  std::size_t diff_order = 1;
  Eigen::VectorXd diff_coeff(precision_order + diff_order);
  Eigen::VectorXd diff_eval_points(4);
  differ_central(1.0e-6, diff_order, precision_order, diff_coeff.data(),
                 diff_eval_points.data());

  for (std::size_t uici = 0; uici < precision_order + diff_order; uici++) {
    x(0) = Ts_center + diff_eval_points(uici);
    jac_nom.col(0) += diff_coeff(uici) * cnstrt.__GetValues(x);
  }

  /*
  x(0) = Ts_center;
  for (std::size_t uici = 0; uici < precision_order + diff_order; uici++) {
    x(1) = sf_center + diff_eval_points(uici);
    jac_nom.col(1) += diff_coeff(uici) * cnstrt.__GetValues(x);
  }
*/

  std::cout << "\n ----- \n" << jac_test << "\n ----- \n";

  Eigen::MatrixXd err_mat = (jac_nom - jac_test).array().abs().matrix();

  assert(err_mat.maxCoeff() < 1.0e-9);
  return 0;
}

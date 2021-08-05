#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines++/ipopt_solver.hpp>
#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>
#include <opstop/torque_constraint.hpp>
#include <pinocchio/parsers/urdf.hpp>

#include <opstop/differ.hpp>

std::size_t number_of_wp = 3;
std::size_t codom_dim = 7;
std::size_t nglp = 10;

double exec_time = (double)number_of_wp - 1.0;
double T = exec_time * 0.9;
double ti = 0.5 * exec_time;
double Ts = 0.8 * T - ti;
double sf = ti * (Ts / ti + 1);

double eta_0 = 1.0 / 15.0 * (exec_time / ti - 1.0);
double xi_0 = (exec_time / ti - 1.0);

double sf_center = ((eta_0 - exec_time / ti) / 2.0 + eta_0) * ti;
double sf_radius = ((eta_0 - exec_time / ti) / 2.0) * ti;
double Ts_center = xi_0 / 2.0 * ti;
double Ts_radius = xi_0 / 2.0 * ti;

Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));

gsplines::PiecewiseFunction trj = gsplines_opt::optimal_sobolev_norm(
    wp, gsplines::basis::BasisLegendre(6), {{4, 1}}, exec_time);

Eigen::VectorXd pol_coeff(6);

int main() {

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
  x(0) = sf_center;
  x(1) = Ts_center;
  Eigen::VectorXd tau = cnstrt.__GetValues(x);

  return 0;
}

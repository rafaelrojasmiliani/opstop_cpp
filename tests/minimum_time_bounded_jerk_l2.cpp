#include <chrono> // for high_resolution_clock
#include <ctime>
#include <fenv.h>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <opstop/ipopt_problem.hpp>
#include <pinocchio/parsers/urdf.hpp>

TEST(Jerk_l2_constraints, stopping) {

  using namespace opstop;
  std::size_t number_of_wp = 10;
  std::size_t codom_dim = 7;
  std::size_t nglp = 10;
  double exec_time = 3 * ((double)number_of_wp - 1.0);
  double ti = 0.5 * exec_time;

  pinocchio::Model model;
  pinocchio::urdf::buildModel("urdf/panda_arm.urdf", model);

  Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));

  gsplines::GSpline trj = gsplines::optimization::optimal_sobolev_norm(
      wp, gsplines::basis::BasisLegendre(6), {{1.0, 3}}, exec_time);
  gsplines::functions::FunctionExpression diffeo = minimum_time_bounded_jerk_l2(
      trj, ti, gsplines::functional_analysis::l2_norm(trj.derivate(3)), model,
      model.effortLimit, 5);
  /*

    gsplines::functions::FunctionExpression diffeo_diff_1 = diffeo.derivate();
    gsplines::functions::FunctionExpression diffeo_diff_2 =
        diffeo_diff_1.derivate();

    double T = diffeo.get_domain().second;

    Eigen::VectorXd time_span = Eigen::VectorXd::LinSpaced(10, ti, T);

    Eigen::MatrixXd diffeo_hist = diffeo(time_span);

    Eigen::MatrixXd diffeo_diff_1_hist = diffeo_diff_1(time_span);
    Eigen::MatrixXd diffeo_diff_2_hist = diffeo_diff_2(time_span);

    Eigen::VectorXd time_spam(1);
    time_spam(0) = ti;

    double ti_test = diffeo(time_spam)(0, 0);
    double ve = diffeo_diff_1(time_spam)(0, 0);
    double ac = diffeo_diff_2(time_spam)(0, 0);

    EXPECT_LT(std::abs(ti_test - ti), 1.0e-9);
    EXPECT_LT(std::abs(ve - 1.0), 1.0e-9);
    EXPECT_LT(std::abs(ac), 1.0e-9);

    gsplines::functions::FunctionExpression stop_trj = trj.compose(diffeo);
    gsplines::functions::FunctionExpression stop_trj_diff_1 =
    stop_trj.derivate(); gsplines::functions::FunctionExpression stop_trj_diff_2
    = stop_trj.derivate(2);

    time_spam(0) = T;
    ti_test = diffeo(time_spam)(0, 0);
    ve = diffeo_diff_1(time_spam)(0, 0);
    ac = diffeo_diff_2(time_spam)(0, 0);

    EXPECT_LT(std::abs(ve), 1.0e-9);
    EXPECT_LT(std::abs(ac), 1.0e-9);

    double vel_value = stop_trj_diff_1(time_spam).row(0).norm();
    double acc_value = stop_trj_diff_2(time_spam).row(0).norm();

    EXPECT_LT(vel_value, 1.0e-9);
    EXPECT_LT(acc_value, 1.0e-9);
    */
}

int main(int argc, char **argv) {
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
  return 0;
}

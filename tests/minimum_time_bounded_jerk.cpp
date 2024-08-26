#include <fenv.h>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gtest/gtest.h>
#include <opstop/ipopt_problem.hpp>

using namespace opstop;

TEST(MINIMUN_TIME_BOUNDED_Jerk, Compuation) {
  std::size_t number_of_wp = 3;
  std::size_t codom_dim = 7;
  std::size_t nglp = 10;

  double exec_time = (double)number_of_wp - 1.0;
  double ti = 0.5 * exec_time;

  Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));

  gsplines::GSpline trj =
      gsplines::optimization::optimal_sobolev_norm(
          wp, gsplines::basis::BasisLegendre(6), {{1.0, 3}}, exec_time)
          .value();

  gsplines::functions::FunctionExpression diffeo =
      minimum_time_bounded_jerk(trj, ti, 100).value();

  gsplines::functions::FunctionExpression diffeo_diff_1 = diffeo.derivate();
  gsplines::functions::FunctionExpression diffeo_diff_2 =
      diffeo_diff_1.derivate();

  double T = diffeo.get_domain().second;

  Eigen::VectorXd time_span = Eigen::VectorXd::LinSpaced(10, ti, T);

  Eigen::MatrixXd diffeo_hist = diffeo(time_span);

  Eigen::MatrixXd diffeo_diff_1_hist = diffeo_diff_1(time_span);
  Eigen::MatrixXd diffeo_diff_2_hist = diffeo_diff_2(time_span);

  Eigen::VectorXd vec(1);
  /// -----------------------------------------
  /// -----------------------------------------
  /// Test of the diffeo
  /// -----------------------------------------
  /// -----------------------------------------
  /// -----------------------------------------
  /// Test stopping parametrization at the beginning
  /// of the stopping procedure
  /// --------------------------------------
  vec(0) = ti;
  ///
  /// Here we check that the diffeo evaluated at ti is equal to ti, this is
  /// because the the diffeo evaluated at ti will give us the values of the
  /// trajectory at ti so q \circ s(t_i) = q(t_i)
  double ti_test = diffeo(vec)(0, 0);
  EXPECT_NEAR(ti_test, ti, 1.0e-7);

  /// Here we check that the derivative of the diffeo is one.
  /// This is because then continuity condition between the original
  /// and the stoping trajectory.
  /// trajectory at ti so {d / d t {q \circ s } }(t_i) = \dot q(t_i)
  double ve = diffeo_diff_1(vec)(0, 0);
  EXPECT_NEAR(ve, 1.0, 1.0e-7);

  /// Here we check that the 2nd derivative of the diffeo is one.
  /// This is because then continuity condition between the original
  /// and the stoping trajectory.
  /// trajectory at ti so {d^2 / d t^2 {q \circ s } }(t_i) = \ddot q(t_i)
  double ac = diffeo_diff_2(vec)(0, 0);
  EXPECT_NEAR(ac, 0.0, 1.0e-7);
  /// -----------------------------------------
  /// Test stopping parametrization at the end
  /// of the stopping procedure
  /// --------------------------------------
  vec(0) = T;
  ti_test = diffeo(vec)(0, 0);
  /// Here we check that the topping condition in velocity
  /// {d / d t s }(T) = 0.0
  ve = diffeo_diff_1(vec)(0, 0);
  EXPECT_NEAR(ve, 0.0, 1.0e-7);
  /// Here we check that the topping condition in acceleration
  ac = diffeo_diff_2(vec)(0, 0);
  EXPECT_NEAR(ac, 0.0, 1.0e-7);

  /// -----------------------------------------
  /// -----------------------------------------
  /// Test of the stopping trajectory
  /// -----------------------------------------
  /// -----------------------------------------
  gsplines::functions::FunctionExpression stop_trj = trj.compose(diffeo);
  gsplines::functions::FunctionExpression stop_trj_diff_1 = stop_trj.derivate();
  gsplines::functions::FunctionExpression stop_trj_diff_2 =
      stop_trj.derivate(2);

  vec(0) = T;
  // Test that the stopping trajectory actually stops
  ve = stop_trj_diff_1(vec).row(0).norm();
  EXPECT_NEAR(ve, 0.0, 1.0e-7);
  // Test that the stopping trajectory actually stops with zero acceleration
  ac = stop_trj_diff_2(vec).row(0).norm();
  EXPECT_NEAR(ac, 0.0, 1.0e-7);
}
int main(int argc, char **argv) {
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

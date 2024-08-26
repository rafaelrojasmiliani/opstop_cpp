#ifndef IPOPT_PROBLEM
#define IPOPT_PROBLEM
#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>
#include <opstop/acceleration_constraints.hpp>
#include <opstop/diffeo_constraints.hpp>
#include <opstop/excursion_cost.hpp>
#include <opstop/jerk_constraints.hpp>
#include <opstop/jerk_l2_constraints.hpp>
#include <opstop/parametrization_variables.hpp>
#include <opstop/time_cost.hpp>
#include <opstop/torque_constraint.hpp>

#include <gsplines/Functions/FunctionBase.hpp>

#include <optional>

namespace opstop {

namespace optimization {
class IpoptSolverOptions {
private:
  static std::optional<IpoptSolverOptions> instance_;
  IpoptSolverOptions();

  std::vector<std::pair<std::string, std::string>> string_options_ = {
      {"linear_solver", "mumps"},
      {"jacobian_approximation", "exact"},
      {"fast_step_computation", "yes"},
      {"derivative_test", "none"},
      {"hessian_approximation", "limited-memory"},
      {"jac_c_constant", "no"},
      {"print_timing_statistics", "yes"},
      {"dependency_detector", "mumps"},
      {"dependency_detection_with_rhs", "no"}};

  std::vector<std::pair<std::string, int>> int_options_ = {{"print_level", 5}};

  std::vector<std::pair<std::string, double>> double_options_ = {
      {"tol", 1.0e-3}};

public:
  static IpoptSolverOptions &instance();

  static void set_option(const std::string &_option_name,
                         const std::string &_option_value);

  static void set_option(const std::string &_option_name, int _option_value);

  static void set_option(const std::string &_option_name, double _option_value);

  static void set_options_on_interface(ifopt::IpoptSolver &solver);
};
} // namespace optimization
class TimeOptimalStopProblem {
private:
  ifopt::Problem nlp_;
  ifopt::IpoptSolver solver_;

public:
  TimeOptimalStopProblem();
};

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_acceleration(const gsplines::functions::FunctionBase &_trj,
                                  double _ti, double _alpha,
                                  const pinocchio::Model &_model,
                                  std::size_t _nglp);

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_jerk_l2(const gsplines::functions::FunctionBase &_trj,
                             double _ti, double _alpha,
                             const pinocchio::Model &_model, std::size_t _nglp);

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_acceleration(const gsplines::functions::FunctionBase &_trj,
                                  double _ti,
                                  const Eigen::VectorXd &_acc_bounds,
                                  const pinocchio::Model &_model,
                                  const Eigen::VectorXd &_torque_bounds,
                                  std::size_t _nglp);

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_jerk(const gsplines::functions::FunctionBase &_trj,
                          double _ti, double _acc_bound);

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_jerk(const gsplines::functions::FunctionBase &_trj,
                          double _ti, std::vector<double> _acc_bounds);

ifopt::Problem
base_minimum_time_problem(const gsplines::functions::FunctionBase &_trj,
                          double _ti);

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_jerk_l2(const gsplines::functions::FunctionBase &_trj,
                             double _ti, double _jerk_l2_bound,
                             pinocchio::Model _model,
                             const Eigen::VectorXd &_torque_bounds,
                             std::size_t _nglp);
} // namespace opstop

#endif /* ifndef IPOPT_PROBLEM                                                 \
        */

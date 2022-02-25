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

namespace opstop {

class TimeOptimalStopProblem {
private:
  ifopt::Problem nlp_;
  ifopt::IpoptSolver solver_;

public:
  TimeOptimalStopProblem();
};

gsplines::functions::FunctionExpression
minimum_time_bouded_acceleration(const gsplines::functions::FunctionBase &_trj,
                                 double _ti, double _acc_bound,
                                 pinocchio::Model _model);

gsplines::functions::FunctionExpression minimum_time_bouded_acceleration(
    const gsplines::functions::FunctionBase &_trj, double _ti,
    const Eigen::VectorXd &_acc_bounds, pinocchio::Model _model,
    const Eigen::VectorXd &_torque_bounds, std::size_t _nglp);

gsplines::functions::FunctionExpression
minimum_time_bouded_jerk(const gsplines::functions::FunctionBase &_trj,
                         double _ti, double _acc_bound);

gsplines::functions::FunctionExpression
minimum_time_bouded_jerk(const gsplines::functions::FunctionBase &_trj,
                         double _ti, std::vector<double> _acc_bounds);

ifopt::Problem
base_minimum_time_problem(const gsplines::functions::FunctionBase &_trj,
                          double _ti);

gsplines::functions::FunctionExpression minimum_time_bouded_jerk_l2(
    const gsplines::functions::FunctionBase &_trj, double _ti,
    double _jerk_l2_bound, pinocchio::Model _model,
    const Eigen::VectorXd &_torque_bounds, std::size_t _nglp);
} // namespace opstop

#endif /* ifndef IPOPT_PROBLEM                                                 \
        */

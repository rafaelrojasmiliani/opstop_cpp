#include <gsplines/FunctionalAnalysis/Sobolev.hpp>
#include <opstop/ipopt_problem.hpp>
#include <sys/resource.h>
#include <unistd.h>
namespace opstop {

ifopt::Problem
base_minimum_time_problem(const gsplines::functions::FunctionBase &_trj,
                          double _ti) {

  std::shared_ptr<ParametrizationVariables> variable =
      std::make_shared<ParametrizationVariables>(_ti, _trj.get_domain().second);

  std::shared_ptr<DiffeoConstraints> diffeo_con =
      std::make_shared<DiffeoConstraints>(_ti, _trj.get_domain().second);

  std::shared_ptr<TimeCost> cost_function = std::make_shared<TimeCost>();

  ifopt::Problem nlp;

  nlp.AddVariableSet(variable);
  nlp.AddConstraintSet(diffeo_con);
  nlp.AddCostSet(cost_function);

  return nlp;
}
/*
gsplines::functions::FunctionExpression
minimum_time_bouded_acceleration(const gsplines::functions::FunctionBase &_trj,
                                 double _ti, double _acc_bound) {
  std::vector<double> bounds(_trj.get_codom_dim(), _acc_bound);
  return minimum_time_bouded_acceleration(_trj, _ti, bounds);
}
*/
gsplines::functions::FunctionExpression minimum_time_bouded_acceleration(
    const gsplines::functions::FunctionBase &_trj, double _ti,
    const Eigen::VectorXd &_acc_bounds, pinocchio::Model _model,
    const Eigen::VectorXd &_torque_bounds, std::size_t _nglp) {

  ifopt::Problem nlp = base_minimum_time_problem(_trj, _ti);

  std::shared_ptr<AccelerationConstraints> acc_con =
      std::make_shared<AccelerationConstraints>(
          _trj, _nglp, _ti,
          std::vector<double>(_acc_bounds.data(),
                              _acc_bounds.data() + _acc_bounds.size()));

  std::shared_ptr<TorqueConstraint> torque_con =
      std::make_shared<TorqueConstraint>(
          _trj, _nglp, _ti,
          std::vector<double>(_torque_bounds.data(),
                              _torque_bounds.data() + _torque_bounds.size()),
          _model);

  // 2. Use the problem objects to build the problem
  nlp.AddConstraintSet(acc_con);
  nlp.AddConstraintSet(torque_con);
  printf("nice = %d\n", getpriority(PRIO_PROCESS, getpid()));
  // nlp.PrintCurrent();
  /*
    printf("-----------------------------\n");
    printf("Number of variables %i \n", nlp.GetNumberOfOptimizationVariables());
    printf("Number of constraints %i \n", nlp.GetNumberOfConstraints());
    printf("-----------------------------\n");
  */
  // 3. Instantiate ipopt solver
  ifopt::IpoptSolver ipopt;
  // 3.1 Customize the solver
  ipopt.SetOption("linear_solver", "ma27");
  ipopt.SetOption("fast_step_computation", "yes");
  ipopt.SetOption("jacobian_approximation", "exact");
  ipopt.SetOption("hessian_approximation", "limited-memory");
  ipopt.SetOption("tol", 1.0e-2);
  ipopt.SetOption("print_level", 5);
  ipopt.SetOption("print_timing_statistics", "yes");
  ipopt.SetOption("mu_target", 500);

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);
  Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();

  double tf = _trj.get_domain().second;
  // printf("Ts = %lf     sf = %lf  ti = %lf\n", x(0), x(1), _ti);
  printf("Ts = %lf     sf = %lf  ti = %lf eta = %lf xi = %lf\n", x(0), x(1),
         _ti, (x(1) - _ti) / (tf - _ti), x(0) / (tf - _ti));

  return get_diffeo(_ti, x(0), x(1));
}

gsplines::functions::FunctionExpression
minimum_time_bouded_jerk(const gsplines::functions::FunctionBase &_trj,
                         double _ti, double _acc_bound) {
  std::vector<double> bounds(_trj.get_codom_dim(), _acc_bound);
  return minimum_time_bouded_jerk(_trj, _ti, bounds);
}
gsplines::functions::FunctionExpression
minimum_time_bouded_jerk(const gsplines::functions::FunctionBase &_trj,
                         double _ti, std::vector<double> _acc_bounds) {

  ifopt::Problem nlp;

  // 1. Build problem objects
  // 1.1 Variables
  std::shared_ptr<ParametrizationVariables> variable =
      std::make_shared<ParametrizationVariables>(_ti, _trj.get_domain().second);
  // 1.2 Constraints
  std::shared_ptr<DiffeoConstraints> diffeo_con =
      std::make_shared<DiffeoConstraints>(_ti, _trj.get_domain().second);

  std::shared_ptr<JerkConstraints> acc_con =
      std::make_shared<JerkConstraints>(_trj, 10, _ti, _acc_bounds);
  // 1.3 Cost Function
  std::shared_ptr<TimeCost> cost_function = std::make_shared<TimeCost>();

  // 2. Use the problem objects to build the problem
  nlp.AddVariableSet(variable);
  nlp.AddConstraintSet(diffeo_con);
  nlp.AddConstraintSet(acc_con);
  nlp.AddCostSet(cost_function);
  nlp.PrintCurrent();

  // 3. Instantiate ipopt solver
  ifopt::IpoptSolver ipopt;
  // 3.1 Customize the solver
  ipopt.SetOption("linear_solver", "mumps");
  ipopt.SetOption("jacobian_approximation", "exact");
  ipopt.SetOption("hessian_approximation", "limited-memory");

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);
  Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();

  printf("Ts = %lf     sf = %lf  ti = %lf\n", x(0), x(1), _ti);

  return get_diffeo(_ti, x(0), x(1));
}

gsplines::functions::FunctionExpression minimum_time_bouded_jerk_l2(
    const gsplines::functions::FunctionBase &_trj, double _ti,
    double _jerk_l2_bound, pinocchio::Model _model,
    const Eigen::VectorXd &_torque_bounds, std::size_t _nglp) {

  ifopt::Problem nlp = base_minimum_time_problem(_trj, _ti);

  double jerk_l2 =
      gsplines::functional_analysis::l2_norm(_trj.derivate(3), 10, 10);

  std::shared_ptr<JerkL2Constraints> acc_con =
      std::make_shared<JerkL2Constraints>(_trj, _nglp, _ti, jerk_l2);

  std::shared_ptr<TorqueConstraint> torque_con =
      std::make_shared<TorqueConstraint>(
          _trj, _nglp, _ti,
          std::vector<double>(_torque_bounds.data(),
                              _torque_bounds.data() + _torque_bounds.size()),
          _model);

  // 2. Use the problem objects to build the problem
  nlp.AddConstraintSet(acc_con);
  nlp.AddConstraintSet(torque_con);
  printf("nice = %d\n", getpriority(PRIO_PROCESS, getpid()));
  // nlp.PrintCurrent();
  /*
    printf("-----------------------------\n");
    printf("Number of variables %i \n", nlp.GetNumberOfOptimizationVariables());
    printf("Number of constraints %i \n", nlp.GetNumberOfConstraints());
    printf("-----------------------------\n");
  */
  // 3. Instantiate ipopt solver
  ifopt::IpoptSolver ipopt;
  // 3.1 Customize the solver
  ipopt.SetOption("linear_solver", "ma27");
  ipopt.SetOption("fast_step_computation", "yes");
  ipopt.SetOption("jacobian_approximation", "exact");
  ipopt.SetOption("hessian_approximation", "limited-memory");
  ipopt.SetOption("tol", 1.0e-2);
  ipopt.SetOption("print_level", 5);
  ipopt.SetOption("print_timing_statistics", "yes");
  ipopt.SetOption("mu_target", 500);

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);
  Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();

  double tf = _trj.get_domain().second;
  // printf("Ts = %lf     sf = %lf  ti = %lf\n", x(0), x(1), _ti);
  printf("Ts = %lf     sf = %lf  ti = %lf eta = %lf xi = %lf\n", x(0), x(1),
         _ti, (x(1) - _ti) / (tf - _ti), x(0) / (tf - _ti));

  return get_diffeo(_ti, x(0), x(1));
}
} // namespace opstop

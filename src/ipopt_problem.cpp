#include <gsplines/FunctionalAnalysis/Sobolev.hpp>
#include <opstop/ipopt_problem.hpp>
#include <sys/resource.h>
#include <unistd.h>
namespace opstop {

namespace optimization {
enum class IpoptReturnStatus : int {
  Solve_Succeeded = 0,
  Solved_To_Acceptable_Level = 1,
  Infeasible_Problem_Detected = 2,
  Search_Direction_Becomes_Too_Small = 3,
  Diverging_Iterates = 4,
  User_Requested_Stop = 5,
  Feasible_Point_Found = 6,

  Maximum_Iterations_Exceeded = -1,
  Restoration_Failed = -2,
  Error_In_Step_Computation = -3,
  Maximum_CpuTime_Exceeded = -4,
  Not_Enough_Degrees_Of_Freedom = -10,
  Invalid_Problem_Definition = -11,
  Invalid_Option = -12,
  Invalid_Number_Detected = -13,

  Unrecoverable_Exception = -100,
  NonIpopt_Exception_Thrown = -101,
  Insufficient_Memory = -102,
  Internal_Error = -199
};

std::optional<IpoptSolverOptions> IpoptSolverOptions::instance_ = std::nullopt;

IpoptSolverOptions::IpoptSolverOptions() = default;

IpoptSolverOptions &IpoptSolverOptions::instance() {
  if (!instance_.has_value()) {
    instance_ = IpoptSolverOptions();
  }
  return instance_.value();
}

void IpoptSolverOptions::set_option(const std::string &_option_name,
                                    const std::string &_option_value) {
  auto iter = std::find_if(
      instance().string_options_.begin(), instance().string_options_.end(),
      [_option_name](const auto &in) { return in.first == _option_name; });

  if (iter == instance().string_options_.end()) {
    instance().string_options_.emplace_back(_option_name, _option_value);
  } else {
    iter->second = _option_value;
  }
}

void IpoptSolverOptions::set_option(const std::string &_option_name,
                                    int _option_value) {
  auto iter = std::find_if(
      instance().int_options_.begin(), instance().int_options_.end(),
      [_option_name](const auto &in) { return in.first == _option_name; });

  if (iter == instance().int_options_.end()) {
    instance().int_options_.emplace_back(_option_name, _option_value);
  } else {
    iter->second = _option_value;
  }
}
void IpoptSolverOptions::set_option(const std::string &_option_name,
                                    double _option_value) {
  auto iter = std::find_if(
      instance().double_options_.begin(), instance().double_options_.end(),
      [_option_name](const auto &in) { return in.first == _option_name; });

  if (iter == instance().double_options_.end()) {
    instance().double_options_.emplace_back(_option_name, _option_value);
  } else {
    iter->second = _option_value;
  }
}
void IpoptSolverOptions::set_options_on_interface(ifopt::IpoptSolver &solver) {
  for (const auto &p : instance().string_options_) {
    solver.SetOption(p.first, p.second);
  }

  for (const auto &p : instance().int_options_) {
    solver.SetOption(p.first, p.second);
  }
  for (const auto &p : instance().double_options_) {
    solver.SetOption(p.first, p.second);
  }
}
} // namespace optimization
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

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_acceleration(const gsplines::functions::FunctionBase &_trj,
                                  double _ti, double _alpha,
                                  const pinocchio::Model &_model,
                                  std::size_t _nglp) {
  std::size_t number_of_segments = 100;
  Eigen::VectorXd bounds =
      _alpha * _trj.derivate(2)
                   .value(Eigen::VectorXd::LinSpaced(number_of_segments + 1,
                                                     _trj.get_domain().first,
                                                     _trj.get_domain().second))
                   .array()
                   .abs()
                   .colwise()
                   .maxCoeff();
  std::cout << "acceleration bound :\n" << bounds << "\n ...\n";
  return minimum_time_bounded_acceleration(_trj, _ti, bounds, _model,
                                           _model.effortLimit, _nglp);
}

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_jerk_l2(const gsplines::functions::FunctionBase &_trj,
                             double _ti, double _alpha,
                             const pinocchio::Model &_model,
                             std::size_t _nglp) {
  std::size_t number_of_segments = 100;
  double jerk_bound =
      _alpha *
      (gsplines::functional_analysis::l2_norm(_trj.derivate(3)) -
       gsplines::functional_analysis::l2_norm(_trj.derivate(3).compose(
           gsplines::functions::Identity({_ti, _trj.get_domain().second}))));
  return minimum_time_bounded_jerk_l2(_trj, _ti, jerk_bound, _model,
                                      _model.effortLimit, _nglp);
}

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_acceleration(const gsplines::functions::FunctionBase &_trj,
                                  double _ti,
                                  const Eigen::VectorXd &_acc_bounds,
                                  const pinocchio::Model &_model,
                                  const Eigen::VectorXd &_torque_bounds,
                                  std::size_t _nglp) {

  std::vector<double> bound(_acc_bounds.data(),
                            _acc_bounds.data() + _acc_bounds.size());

  ifopt::Problem nlp = base_minimum_time_problem(_trj, _ti);

  std::shared_ptr<AccelerationConstraints> acc_con =
      std::make_shared<AccelerationConstraints>(_trj, _nglp, _ti, bound);

  std::shared_ptr<TorqueConstraint> torque_con =
      std::make_shared<TorqueConstraint>(
          _trj, _nglp, _ti,
          std::vector<double>(_torque_bounds.data(),
                              _torque_bounds.data() + _torque_bounds.size()),
          _model);

  // 2. Use the problem objects to build the problem
  nlp.AddConstraintSet(acc_con);
  nlp.AddConstraintSet(torque_con);
  nlp.PrintCurrent();
  /*
    printf("-----------------------------\n");
    printf("Number of variables %i \n", nlp.GetNumberOfOptimizationVariables());
    printf("Number of constraints %i \n", nlp.GetNumberOfConstraints());
    printf("-----------------------------\n");
  */
  // 3. Instantiate ipopt solver
  ifopt::IpoptSolver ipopt;

  // 3.1 Customize the solver

  optimization::IpoptSolverOptions::instance().set_options_on_interface(ipopt);

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);
  if (ipopt.GetReturnStatus() !=
      static_cast<int>(optimization::IpoptReturnStatus::Solve_Succeeded)) {
    return std::nullopt;
  }
  Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();

  double tf = _trj.get_domain().second;
  // printf("Ts = %lf     sf = %lf  ti = %lf\n", x(0), x(1), _ti);
  printf("Ts = %lf     sf = %lf  ti = %lf eta = %lf xi = %lf\n", x(0), x(1),
         _ti, (x(1) - _ti) / (tf - _ti), x(0) / (tf - _ti));

  optimization::IpoptSolverOptions::Solution sol;
  sol.ti = _ti;
  sol.Ts = x(0);
  sol.sf = x(1);
  optimization::IpoptSolverOptions::instance().last_solution.emplace(sol);

  return get_diffeo(_ti, x(0), x(1));
}

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_jerk(const gsplines::functions::FunctionBase &_trj,
                          double _ti, double _acc_bound) {
  std::vector<double> bounds(_trj.get_codom_dim(), _acc_bound);
  return minimum_time_bounded_jerk(_trj, _ti, bounds);
}
std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_jerk(const gsplines::functions::FunctionBase &_trj,
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
  // nlp.PrintCurrent();

  // 3. Instantiate ipopt solver
  ifopt::IpoptSolver ipopt;
  // 3.1 Customize the solver
  optimization::IpoptSolverOptions::instance().set_options_on_interface(ipopt);

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);
  if (ipopt.GetReturnStatus() !=
      static_cast<int>(optimization::IpoptReturnStatus::Solve_Succeeded)) {
    return std::nullopt;
  }
  Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();

  printf("Ts = %lf     sf = %lf  ti = %lf\n", x(0), x(1), _ti);

  optimization::IpoptSolverOptions::Solution sol;
  sol.ti = _ti;
  sol.Ts = x(0);
  sol.sf = x(1);
  optimization::IpoptSolverOptions::instance().last_solution.emplace(sol);
  return get_diffeo(_ti, x(0), x(1));
}

std::optional<gsplines::functions::FunctionExpression>
minimum_time_bounded_jerk_l2(const gsplines::functions::FunctionBase &_trj,
                             double _ti, double _jerk_l2_bound,
                             pinocchio::Model _model,
                             const Eigen::VectorXd &_torque_bounds,
                             std::size_t _nglp) {

  ifopt::Problem nlp = base_minimum_time_problem(_trj, _ti);

  double jerk_l2 =
      gsplines::functional_analysis::l2_norm(_trj.derivate(3), 10, 10);

  std::shared_ptr<JerkL2Constraints> jerk_l2_con =
      std::make_shared<JerkL2Constraints>(_trj, _nglp, _ti, jerk_l2);

  std::shared_ptr<TorqueConstraint> torque_con =
      std::make_shared<TorqueConstraint>(
          _trj, _nglp, _ti,
          std::vector<double>(_torque_bounds.data(),
                              _torque_bounds.data() + _torque_bounds.size()),
          _model);

  // 2. Use the problem objects to build the problem
  nlp.AddConstraintSet(jerk_l2_con);
  nlp.AddConstraintSet(torque_con);
  // nlp.PrintCurrent();
  // 3. Instantiate ipopt solver
  ifopt::IpoptSolver ipopt;
  // 3.1 Customize the solver
  optimization::IpoptSolverOptions::instance().set_options_on_interface(ipopt);

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);

  if (ipopt.GetReturnStatus() !=
      static_cast<int>(optimization::IpoptReturnStatus::Solve_Succeeded)) {
    return std::nullopt;
  }
  Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();

  double tf = _trj.get_domain().second;
  // printf("Ts = %lf     sf = %lf  ti = %lf\n", x(0), x(1), _ti);
  printf("Ts = %lf     sf = %lf  ti = %lf eta = %lf xi = %lf\n", x(0), x(1),
         _ti, (x(1) - _ti) / (tf - _ti), x(0) / (tf - _ti));

  return get_diffeo(_ti, x(0), x(1));
}
} // namespace opstop

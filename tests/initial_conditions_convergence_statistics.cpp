#include <Eigen/StdVector>
#include <chrono> // for high_resolution_clock
#include <ctime>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <iostream>
#include <opstop/ipopt_problem.hpp>
#include <pinocchio/parsers/urdf.hpp>
#include <random>

using namespace opstop;

std::random_device rd;
std::mt19937 gen(rd());

std::uniform_real_distribution<float> joint_dis(-3.14, 3.14);

bool test_constraints(ifopt::Problem &_problem) {
  ifopt::Component::VecBound bounds_vector(_problem.GetBoundsOnConstraints());
  Eigen::VectorXd cnstr(
      _problem.EvaluateConstraints(_problem.GetVariableValues().data()));

  for (std::size_t i = 0; i < _problem.GetNumberOfConstraints(); i++) {
    if (cnstr(i) < bounds_vector[i].lower_ or
        cnstr(i) > bounds_vector[i].upper_)
      return false;
  }
  return true;
}

std::tuple<Eigen::VectorXd, int, int, std::vector<Eigen::VectorXd>>
solve_problem(const gsplines::functions::FunctionBase &_trj, double _ti,
              std::size_t _nglp, Eigen::VectorXd &_initial_guess);

bool almost_equal(Eigen::VectorXd &_m_nom, Eigen::VectorXd &_m_test,
                  double _tol) {

  if (_m_nom.array().abs().maxCoeff() < _tol) {
    return ((_m_nom - _m_test).array().abs().maxCoeff() < _tol);
  } else {
    double err = (_m_nom - _m_test).array().abs().maxCoeff() /
                 _m_nom.array().abs().maxCoeff();

    if (err > _tol) {
      return false;
    }
  }
  return true;
}

gsplines::GSpline get_random_curve(std::size_t _codon_dim) {

  std::size_t number_of_wp = 5;
  std::size_t codom_dim = _codon_dim;

  double exec_time = ((double)number_of_wp - 1.0);

  Eigen::MatrixXd wp = Eigen::MatrixXd::NullaryExpr(
      number_of_wp, codom_dim, [&]() { return joint_dis(gen); });

  return gsplines::optimization::optimal_sobolev_norm(
      wp, gsplines::basis::BasisLegendre(6), {{1.0, 3}}, exec_time);
}

Eigen::VectorXd max_acc(const gsplines::functions::FunctionBase &_trj) {
  std::unique_ptr<gsplines::functions::FunctionBase> acc = _trj.deriv(2);

  Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(
      100, _trj.get_domain().first, _trj.get_domain().second);

  Eigen::MatrixXd acc_val = acc->value(time_spam);

  Eigen::VectorXd result(acc_val.cols());

  for (std::size_t i = 0; i < acc_val.cols(); i++) {
    result(i) = acc_val.col(i).array().abs().maxCoeff();
  }

  return result;
}

struct InitialConditionData {
  int number_of_attempts = 0;
  int number_of_successes = 0;
  int max_number_of_iterations = 0;
};

struct ProblemData {
  typedef std::vector<Eigen::VectorXd,
                      Eigen::aligned_allocator<Eigen::VectorXd>>
      VectorXdVector;
  std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>
      solutions;
  std::vector<std::vector<Eigen::VectorXd>> iteratioins_history;
  std::vector<std::size_t> sic_idx;
  double solution_max_error = 0.0;
};

pinocchio::Model g_model;
int main() {

  pinocchio::urdf::buildModel("urdf/panda_arm.urdf", g_model);

  std::size_t n_test = 30;
  int maximum_number_of_iterations = 0;
  int minimum_number_of_iterations = 100;
  std::size_t nglp = 5;

  std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>
      initial_guesses;

  for (double xi = 0.01; xi <= 1.0; xi += 0.1) {
    for (double eta : {2.0 / 5.0 * xi, 3.0 / 5.0 * xi, 1.0 / 2.0 * xi}) {
      Eigen::VectorXd initial_guess(2);
      initial_guess(0) = xi;
      initial_guess(1) = eta;
      initial_guesses.push_back(initial_guess);
    }
  }

  std::vector<std::tuple<gsplines::GSpline, double, ProblemData>>
      random_problems;

  for (std::size_t i = 0; i < n_test; i++) {
    gsplines::GSpline trj = get_random_curve(7);

    Eigen::VectorXd acc_bounds = 5.0 * max_acc(trj);
    std::uniform_real_distribution<float> ti_dis(trj.get_domain().first,
                                                 trj.get_domain().second);

    for (std::size_t j = 0; j < 0.3 * n_test; j++) {
      double ti = ti_dis(gen);
      random_problems.push_back(std::make_tuple(trj, ti, ProblemData()));
    }
  }

  Eigen::VectorXd solution;

  std::size_t pc = 0; // problem counter
  for (auto &problem : random_problems) {
    std::size_t igc = 0; // initial guess counter
    ProblemData &problem_data = std::get<2>(problem);
    for (Eigen::VectorXd &initial_guess : initial_guesses) {
      int status, iterations;
      problem_data.iteratioins_history.emplace_back(0);
      std::tie(solution, status, iterations,
               problem_data.iteratioins_history.back()) =
          solve_problem(std::get<0>(problem), std::get<1>(problem), 5,
                        initial_guess);

      if (iterations > -1 and status == 0) {
        problem_data.sic_idx.push_back(igc);

        for (auto &x : problem_data.solutions) {
          double err = (x - solution).array().abs().maxCoeff();
          if (err > problem_data.solution_max_error) {
            problem_data.solution_max_error = err;
          }
        }
        problem_data.solutions.push_back(solution);
      } else {
        problem_data.iteratioins_history.pop_back();
      }
      igc++;
    }
    pc++; // problem counter
  }

  for (auto &problem : random_problems) {
    std::size_t igc = 0; // initial guess counter
    ProblemData &problem_data = std::get<2>(problem);
    printf("maximum error %14.7e \n", problem_data.solution_max_error);
  }

  char file_name[200];
  for (std::size_t i : {0, 100, 200, 400, 20}) {
    std::size_t i2;
    for (i2 = i; i2 < random_problems.size(); i2++) {
      ProblemData &problem_data = std::get<2>(random_problems[i2]);
      if (not problem_data.iteratioins_history.empty())
        break;
    }
    if (i2 == random_problems.size())
      break;
    ProblemData &problem_data = std::get<2>(random_problems[i2]);
    std::size_t iter1 = 0;
    for (std::size_t ic = 0; ic < problem_data.iteratioins_history.size();
         ic++) {
      sprintf(file_name, "iterations_history_%02zu_%02zu_ti=%.4lf.txt", i2, ic,
              std::get<1>(random_problems[i2]));
      FILE *output = fopen(file_name, "w");
      for (Eigen::VectorXd &vec : problem_data.iteratioins_history[ic]) {
        fprintf(output, "%+14.7e %+14.7e \n", vec(0), vec(1));
      }
      fclose(output);
    }
  }
  return 0;
}

std::tuple<Eigen::VectorXd, int, int, std::vector<Eigen::VectorXd>>
solve_problem(const gsplines::functions::FunctionBase &_trj, double _ti,
              std::size_t _nglp, Eigen::VectorXd &_initial_guess) {

  std::tuple<Eigen::VectorXd, int, int> result = {Eigen::VectorXd(2), -1, -1};
  Eigen::VectorXd solution(2);

  // Initial guess correction
  Eigen::VectorXd initial_guess(2);
  initial_guess(0) = _initial_guess(0) * (_trj.get_domain_length() - _ti);
  initial_guess(1) = _initial_guess(1) * (_trj.get_domain_length() - _ti) + _ti;
  // ifopt problem construction
  std::shared_ptr<ParametrizationVariables> variable =
      std::make_shared<ParametrizationVariables>(_ti, _trj.get_domain().second);

  std::shared_ptr<DiffeoConstraints> diffeo_con =
      std::make_shared<DiffeoConstraints>(_ti, _trj.get_domain().second);

  std::shared_ptr<TimeCost> cost_function = std::make_shared<TimeCost>();

  Eigen::VectorXd acc_bounds = 5.0 * max_acc(_trj);
  std::shared_ptr<AccelerationConstraints> acc_con =
      std::make_shared<AccelerationConstraints>(
          _trj, _nglp, _ti,
          std::vector<double>(acc_bounds.data(),
                              acc_bounds.data() + acc_bounds.size()));

  std::shared_ptr<TorqueConstraint> torque_con =
      std::make_shared<TorqueConstraint>(
          _trj, _nglp, _ti,
          std::vector<double>(g_model.effortLimit.data(),
                              g_model.effortLimit.data() +
                                  g_model.effortLimit.size()),
          g_model);

  variable->SetVariables(initial_guess);
  ifopt::Problem nlp;
  nlp.AddVariableSet(variable);
  nlp.AddConstraintSet(acc_con);
  nlp.AddConstraintSet(diffeo_con);
  nlp.AddCostSet(cost_function);

  if (test_constraints(nlp)) {
    ifopt::IpoptSolver ipopt;
    ipopt.SetOption("linear_solver", "ma27");
    ipopt.SetOption("fast_step_computation", "yes");
    ipopt.SetOption("jacobian_approximation", "exact");
    ipopt.SetOption("hessian_approximation", "limited-memory");
    ipopt.SetOption("tol", 1.0e-5);
    ipopt.SetOption("print_level", 0);
    ipopt.SetOption("max_iter", 100);
    ipopt.Solve(nlp);

    solution = nlp.GetOptVariables()->GetValues();

    return std::make_tuple(solution, ipopt.GetReturnStatus(),
                           nlp.GetIterationCount(), nlp.GetIterations());
  }
  return std::make_tuple(solution, -1, -1, std::vector<Eigen::VectorXd>(0));
}

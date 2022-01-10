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
  int number_of_attractors = 1;
  int number_of_successes = 0;
};

int main() {

  pinocchio::Model model;
  pinocchio::urdf::buildModel("urdf/panda_arm.urdf", model);

  std::size_t n_test = 30;
  int maximum_number_of_iterations = 0;
  int minimum_number_of_iterations = 100;
  std::size_t nglp = 5;

  std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>
      initial_guess_list_dimensionless;
  std::vector<InitialConditionData> initial_condition_data_list;

  for (double xi = 0.01; xi <= 1.0; xi += 0.1) {
    for (double eta : {2.0 / 5.0 * xi, 3.0 / 5.0 * xi, 1.0 / 2.0 * xi}) {
      Eigen::VectorXd initial_guess(2);
      initial_guess(0) = xi;
      initial_guess(1) = eta;
      initial_guess_list_dimensionless.push_back(initial_guess);
      initial_condition_data_list.emplace_back();
    }
  }

  for (std::size_t i = 0; i < n_test; i++) {
    gsplines::GSpline trj = get_random_curve(7);

    Eigen::VectorXd acc_bounds = 5.0 * max_acc(trj);
    std::uniform_real_distribution<float> ti_dis(trj.get_domain().first,
                                                 trj.get_domain().second);

    for (std::size_t j = 0; j < 0.3 * n_test; j++) {
      double ti = ti_dis(gen);
      std::shared_ptr<ParametrizationVariables> variable =
          std::make_shared<ParametrizationVariables>(ti,
                                                     trj.get_domain().second);

      std::shared_ptr<DiffeoConstraints> diffeo_con =
          std::make_shared<DiffeoConstraints>(ti, trj.get_domain().second);

      std::shared_ptr<TimeCost> cost_function = std::make_shared<TimeCost>();

      std::shared_ptr<AccelerationConstraints> acc_con =
          std::make_shared<AccelerationConstraints>(
              trj, nglp, ti,
              std::vector<double>(acc_bounds.data(),
                                  acc_bounds.data() + acc_bounds.size()));

      std::shared_ptr<TorqueConstraint> torque_con =
          std::make_shared<TorqueConstraint>(
              trj, nglp, ti,
              std::vector<double>(model.effortLimit.data(),
                                  model.effortLimit.data() +
                                      model.effortLimit.size()),
              model);

      std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>>
          attractors_list;

      for (int initial_guess_counter = 0;
           initial_guess_counter < initial_guess_list_dimensionless.size();
           initial_guess_counter++) {

        Eigen::VectorXd initial_guess_dimensionless(
            initial_guess_list_dimensionless[initial_guess_counter]);
        double T_s =
            initial_guess_dimensionless(0) * (trj.get_exec_time() - ti);
        double s_f =
            initial_guess_dimensionless(1) * (trj.get_exec_time() - ti) + ti;
        Eigen::VectorXd initial_guess(2);
        initial_guess(0) = T_s;
        initial_guess(1) = s_f;
        variable->SetVariables(initial_guess);
        ifopt::Problem nlp;
        nlp.AddVariableSet(variable);
        nlp.AddConstraintSet(acc_con);
        nlp.AddConstraintSet(diffeo_con);
        nlp.AddCostSet(cost_function);

        if (test_constraints(nlp)) {
          initial_condition_data_list[initial_guess_counter]
              .number_of_attempts++;
          ifopt::IpoptSolver ipopt;
          ipopt.SetOption("linear_solver", "ma27");
          ipopt.SetOption("fast_step_computation", "yes");
          ipopt.SetOption("jacobian_approximation", "exact");
          ipopt.SetOption("hessian_approximation", "limited-memory");
          ipopt.SetOption("tol", 1.0e-2);
          ipopt.SetOption("print_level", 0);
          ipopt.SetOption("max_iter", 100);
          ipopt.Solve(nlp);
          Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();

          if (ipopt.GetReturnStatus() == 0) {

            initial_condition_data_list[initial_guess_counter]
                .number_of_successes++;

            if (std::find_if(attractors_list.begin(), attractors_list.end(),
                             [&x](Eigen::VectorXd &_element) {
                               return almost_equal(_element, x, 5.0e-2);
                             }) == attractors_list.end()) {
              attractors_list.push_back(x);
            }
          }
        }
        if (attractors_list.size() >
            initial_condition_data_list[initial_guess_counter]
                .number_of_attractors)
          initial_condition_data_list[initial_guess_counter]
              .number_of_attractors++;
      }

      std::cout << "\n-----------------\n";
      for (auto vec : attractors_list)
        std::cout << vec << "\n";
      std::cout << "\n-----------------\n";
    }
  }

  for (InitialConditionData &data : initial_condition_data_list) {
    printf(" %d  %d %d \n", data.number_of_attempts, data.number_of_attractors,
           data.number_of_successes);
  }
  return 0;
}

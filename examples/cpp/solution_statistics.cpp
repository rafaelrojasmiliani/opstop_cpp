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

std::uniform_real_distribution<float> joint_dis(-1.0, 1.0);

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
int main() {

  pinocchio::Model model;
  pinocchio::urdf::buildModel("urdf/panda_arm.urdf", model);

  std::size_t n_test = 100;
  std::size_t iters = 0;
  double total_time = 0;
  double total_interations = 0;
  double maximum_time = 0;
  int maximum_number_of_iterations = 0;
  int minimum_number_of_iterations = 100;
  std::size_t nglp = 5;
  std::size_t time_histogram[11] = {0};
  std::size_t iteration_histogram[11] = {0};
  double time_histogram_classes[10] = {2.5,  5.0,  7.5,  10.0, 12.5,
                                       15.0, 17.5, 20.0, 22.5, 25.0};
  std::size_t iteration_histogram_classes[10] = {5,  10, 15, 20, 25,
                                                 30, 35, 40, 45, 50};

  FILE *output = fopen("solution_statistics.txt", "w");

  fprintf(output, "# tf - ti      |       xi     |     eta      | time [ms]  | "
                  "          Status   | iterations\n");

  for (std::size_t i = 0; i < n_test; i++) {
    gsplines::GSpline trj = get_random_curve(7);

    Eigen::VectorXd acc_bounds = 5.0 * max_acc(trj);
    std::uniform_real_distribution<float> ti_dis(trj.get_domain().first,
                                                 trj.get_domain().second);

    for (std::size_t j = 0; j < 0.3 * n_test; j++) {
      auto t1 = std::chrono::high_resolution_clock::now();
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

      if (variable->set_for_acceleration_bounds(trj, ti, acc_bounds)) {

        ifopt::Problem nlp;
        nlp.AddVariableSet(variable);
        nlp.AddConstraintSet(acc_con);
        nlp.AddConstraintSet(diffeo_con);
        nlp.AddCostSet(cost_function);

        ifopt::IpoptSolver ipopt;
        ipopt.SetOption("linear_solver", "mumps");
        ipopt.SetOption("fast_step_computation", "yes");
        ipopt.SetOption("jacobian_approximation", "exact");
        ipopt.SetOption("hessian_approximation", "limited-memory");
        ipopt.SetOption("tol", 1.0e-2);
        ipopt.SetOption("print_level", 0);
        ipopt.SetOption("max_iter", 100);

        ipopt.Solve(nlp);
        Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> ms_double = t2 - t1;
        double dt = trj.get_domain().second - ti;
        double xi = x(0) / dt;
        double eta = (x(1) - ti) / dt;
        total_time += ms_double.count();
        if (ms_double.count() > maximum_time)
          maximum_time = ms_double.count();
        total_interations += nlp.GetIterationCount();
        if (nlp.GetIterationCount() > maximum_number_of_iterations)
          maximum_number_of_iterations = nlp.GetIterationCount();
        if (nlp.GetIterationCount() < minimum_number_of_iterations)
          minimum_number_of_iterations = nlp.GetIterationCount();

        for (std::size_t uichist = 0; uichist < 10; uichist++) {
          if (ms_double.count() < time_histogram_classes[uichist]) {
            time_histogram[uichist]++;
            break;
          }
        }
        for (std::size_t uichist = 0; uichist < 10; uichist++) {
          if (nlp.GetIterationCount() < iteration_histogram_classes[uichist]) {
            iteration_histogram[uichist]++;
            break;
          }
        }

        iters++;
        fprintf(output, " %+14.7E %+14.7E %+14.7E %+14.7E %14d %14d\n", dt, xi,
                eta, ms_double.count(), ipopt.GetReturnStatus(),
                nlp.GetIterationCount());
        fflush(output);
      }
    }
  }
  fclose(output);
  printf("mean optimization time       = %+14.7e ms\n"
         "maximum time                 = %+14.7e ms\n"
         "mean number of iterations    = %+14.7e \n"
         "maximum number of iterations = %14d \n"
         "minumum number of iterations = %14d \n",
         total_time / iters, maximum_time, total_interations / iters,
         maximum_number_of_iterations, minimum_number_of_iterations);

  output = fopen("time_histogram.txt", "w");
  printf("\n----- Time  Histogram -----\n");
  for (std::size_t uichist = 0; uichist < 10; uichist++) {
    printf("class %4.1lf ms          %6.3lf \n",
           time_histogram_classes[uichist],
           static_cast<double>(time_histogram[uichist]) /
               static_cast<double>(iters));
    fprintf(output, "%4.1lf %14.7e \n", time_histogram_classes[uichist],
            static_cast<double>(time_histogram[uichist]) /
                static_cast<double>(iters));
    fflush(output);
  }
  fclose(output);
  output = fopen("iteration_histogram.txt", "w");
  printf("\n----- Number of Iterations Histogram -----\n");
  for (std::size_t uichist = 0; uichist < 10; uichist++) {
    printf("class %2zu iterations    %6.3lf \n",
           iteration_histogram_classes[uichist],
           static_cast<double>(iteration_histogram[uichist]) / iters);
    fprintf(output, "%2zu %14.7e \n", iteration_histogram_classes[uichist],
            static_cast<double>(iteration_histogram[uichist]) / iters);
    fflush(output);
  }
  fclose(output);
  return 0;
}

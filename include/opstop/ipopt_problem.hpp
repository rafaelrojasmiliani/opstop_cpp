#ifndef IPOPT_PROBLEM
#define IPOPT_PROBLEM
#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>
#include <opstop/acceleration_constraints.hpp>
#include <opstop/diffeo_constraints.hpp>
#include <opstop/parametrization_variables.hpp>
#include <opstop/time_cost.hpp>

#include <gsplines++/Functions/FunctionExpression.hpp>

gsplines::functions::FunctionExpression
minimum_time_bouded_acceleration(gsplines::functions::FunctionExpression &_trj,
                                 double _ti, double _acc_bound);
gsplines::functions::FunctionExpression
minimum_time_bouded_acceleration(gsplines::functions::FunctionExpression &_trj,
                                 double _ti, std::vector<double> _acc_bounds);
#endif /* ifndef IPOPT_PROBLEM                                                 \
        */

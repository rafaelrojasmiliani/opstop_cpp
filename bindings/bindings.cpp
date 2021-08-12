#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <gsplines++/Functions/FunctionExpression.hpp>
#include <opstop/ipopt_problem.hpp>
#include <opstop/parametrization.hpp>

PYBIND11_MODULE(opstop, m) {

  m.def("get_diffeo", &opstop::get_diffeo);
  m.def("minimum_time_bouded_acceleration",
        [](gsplines::functions::FunctionExpression &_trj, double _ti,
           std::vector<double> _acc_bounds) {
          return opstop::minimum_time_bouded_acceleration(_trj, _ti,
                                                          _acc_bounds);
        });
}

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <pinocchio/parsers/urdf.hpp>

#include <gsplines/Functions/FunctionExpression.hpp>
#include <opstop/ipopt_problem.hpp>
#include <opstop/parametrization.hpp>

PYBIND11_MODULE(pyopstop, m) {

  pybind11::module::import("gsplines");
  m.def("get_diffeo", &opstop::get_diffeo);
  m.def("get_diffeo_wrt_tau", &opstop::get_diffeo_wrt_tau);
  m.def("get_tau", &opstop::get_tau);
  m.def("get_tau_inv", &opstop::get_tau_inv);
  m.def("minimum_time_bouded_acceleration",
        [](const gsplines::functions::FunctionBase &_trj, double _ti,
           double _alpha, const std::string &_model_path, std::size_t _nglp) {
          pinocchio::Model model;
          pinocchio::urdf::buildModel(_model_path, model);
          return opstop::minimum_time_bouded_acceleration(_trj, _ti, _alpha,
                                                          model, _nglp);
        });
}

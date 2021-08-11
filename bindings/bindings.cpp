#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include <opstop/parametrization.hpp>

PYBIND11_MODULE(opstop, m) { m.def("get_diffeo", &opstop::get_diffeo); }

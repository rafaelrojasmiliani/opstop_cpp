#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(opstop, m) {

  m.def("optimal_sobolev_norm", &gsplines_opt::optimal_sobolev_norm);
  // Operations
  // m.def("optimal_sobolev_norm", &gsplines_opt::optimal_sobolev_norm);
}

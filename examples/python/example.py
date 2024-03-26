#!/usr/bin/python3
"""
Example of how to get a curve that stops in minimum time with acceleration
bounds.
"""
import gsplines
import gsplines.plot as gplot
import pathlib
import sys
import unittest
import numpy as np

try:
    import opstop
except ImportError:
    import pyopstop as opstop

from gsplines.optimization import minimum_jerk_path


def main():
    """
    Main example
    """

    model_file = pathlib.Path(__file__).absolute(
    ).parents[2].joinpath('tests', 'urdf', 'panda_arm.urdf')

    print(str(model_file))

    dim = 7
    intervals = 4
    waypoints = np.random.rand(intervals+1, dim)
    trj = minimum_jerk_path(waypoints).linear_scaling_new_execution_time(5.0)

    stop_ti = trj.get_domain()[1]*0.6

    optimal_parametrization = opstop.minimum_time_bounded_acceleration(
        trj, stop_ti, 1.5, str(model_file), 8)

    stop_trj = trj.compose(optimal_parametrization)

    gplot.plot_compare([stop_trj, trj], ['red', 'blue'], [
                       'Emergency Stop Trajectory',
                       'Original Trajectory'], _show=True, _up_to_deriv=2)

    # new_exec_time = stop_trj.get_domain()[1]

    # stop_trj_d1 = stop_trj.deriv()
    # stop_trj_d2 = stop_trj.deriv(2)
    # self.assertLess(np.max(np.abs(stop_trj_d1([new_exec_time]))), 1.0e-9)
    # self.assertLess(np.max(np.abs(stop_trj_d2([new_exec_time]))), 1.0e-9)

    # stop_time_span = np.arange(stop_ti, new_exec_time, 0.001)
    # self.assertLess(
    #     np.max(np.abs(stop_trj_d2(stop_time_span))), 6.0*1.1)


if __name__ == '__main__':
    main()

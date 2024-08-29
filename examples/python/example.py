#!/usr/bin/python3
"""
Example of how to get a curve that stops in minimum time with acceleration
bounds.
"""
import gsplines
import gsplines.plot as gplot
import unittest
import numpy as np

try:
    import opstop
except ImportError:
    import pathlib

    cwd = pathlib.Path(__file__).parent.absolute()
    setup_script = pathlib.Path(
        cwd, '../../', 'setup_python_env.py')
    exec(setup_script.read_text())
    try:
        import opstop
    except ImportError:
        import sys
        sys.exit(1)


def main():
    """
    Main example
    """

    model_file = pathlib.Path(__file__).absolute(
    ).parents[2].joinpath('tests', 'urdf', 'panda_arm.urdf')

    waypoints = \
        np.array([[-0.237,  0.034,  0.246,  -1.631, -0.007, 1.705, 0.898],
                  [-0.270, -0.010, 0.246,  -1.631, -0.007, 1.705, 0.898],
                  [-0.286, 0.008,  0.222,  -1.616, 0.024,  1.692, 0.856],
                  [-0.377, 0.123,  0.074,  -1.526, 0.210,  1.614, 0.606],
                  [-0.387, 0.145,  0.058,  -1.517, 0.250,  1.592, 0.565],
                  [-0.491, 0.365,  -0.102, -1.426, 0.646,  1.365, 0.158]]
                 )
    numberOfWaypoints, dimensionOfAmbientSpace = waypoints.shape
    ni = numberOfWaypoints - 1
    execution_time = 4 * ni / np.sqrt(2)
    k = 1.5
    alpha = np.power(k, 4) / (1.0 + np.power(k, 4))
    basis = gsplines.basis.Basis0101(alpha)
    trj = gsplines.optimization.optimal_sobolev_norm(
        waypoints, basis, [(1, alpha), (3, 1-alpha)], execution_time)

    trj = trj.linear_scaling_new_execution_time(2.00053)

    stop_ti = 1.45928

    alpha_acc_bound = 1.98

    optimal_parametrization = opstop.minimum_time_bounded_acceleration(
        trj, stop_ti, alpha_acc_bound, str(model_file), 8)

    stop_trj = trj.compose(optimal_parametrization)

    print(alpha)
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

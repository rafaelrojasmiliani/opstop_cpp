#!/usr/bin/python3
"""
Example of how to get a curve that stops in minimum time with acceleration
bounds.
"""
import gsplines
import gsplines.plot as gplot
import unittest
import numpy as np
import pathlib

try:
    import opstop
except ImportError:

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

# gspline:
#   domain_left_boundary: 0.0
#   domain_right_boundary: 2.164836572
#   codom_dim: 9
#   number_of_intervals: 1
#   basis:
#     name: "legendre"
#     dim: 4
#     parameters: []
#   coefficients: [1.7298899103880383, 0.03457196628719199, -0.0, -0.005761994381198665, 0.22234435729877736, -0.054186120424065586, -0.0, 0.009031020070677598, -1.6216245522995387, -1.0193175307026285, -4.005752344641659e-17, 0.16988625511710476, -2.2492528887271037, -0.2705483863218848, -0.0, 0.04509139772031414, 0.2046680011084518, -0.029006004707747935, -1.2517976077005184e-18, 0.004834334117957989, 2.2778514981148423, 0.06212618368823595, -0.0, -0.010354363948039324, 0.7498029371412586, -0.9926840055358922, -0.0, 0.16544733425598204, -0.0027804829725030394, 0.0, -0.0, 0.0, -0.0027804829725030394, 0.0, -0.0, 0.0]
#   interval_lengths: [2.164836572]

    basis = gsplines.basis.get_basis("legendre",
                                     4, [])

    domain = (0.0,
              2.164836572)

    codom_dim = 9

    number_of_intervals = 1

    coefficients = [0.024982879472156383, -0.26031064300024337, -0.0, 0.043385107166707224, -0.07528467536599337, 0.012468025015078571, 5.111428525809518e-19, -0.0020780041691797614, -0.0025524799879010907, -0.5571324525407222, -0.0, 0.09285540875678704, -1.6204058549081795, 0.002092858579872381, -0.0, -0.00034880976331206346,
                    0.022289364507168237, -0.037609615284698694, -0.0, 0.006268269214116449, 1.567069257430161, -0.0025737185527551, -0.0, 0.0004289530921258499, 0.7891473199553142, -0.8157687812572059, -3.2713142565180915e-17, 0.13596146354286762, -0.002763126353232897, 0.0, -0.0, 0.0, -0.002763126353232897, 0.0, -0.0, 0.0]

    interval_lengths = [1.76792276]

    trj = gsplines.GSpline(domain, codom_dim, number_of_intervals,
                           basis, coefficients, interval_lengths, "gspline")

    model_file = pathlib.Path(__file__).absolute(
    ).parents[2].joinpath('tests', 'urdf', 'panda_arm_with_fingers.urdf')

    # waypoints = \
    #     np.array([[-0.237,  0.034,  0.246,  -1.631, -0.007, 1.705, 0.898, -0.002, -0.002],
    #               [-0.270, -0.010, 0.246,  -1.631, -
    #                   0.007, 1.705, 0.898, -0.002, -0.002],
    #               [-0.286, 0.008,  0.222,  -1.616, 0.024,
    #                   1.692, 0.856, -0.002, -0.002],
    #               [-0.377, 0.123,  0.074,  -1.526, 0.210,
    #                   1.614, 0.606, -0.002, -0.002],
    #               [-0.387, 0.145,  0.058,  -1.517, 0.250,
    #                   1.592, 0.565, -0.002, -0.002],
    #               [-0.491, 0.365,  -0.102, -1.426, 0.646,  1.365, 0.158, -0.002, -0.002]]
    #              )
    # numberOfWaypoints, dimensionOfAmbientSpace = waypoints.shape
    # ni = numberOfWaypoints - 1
    # execution_time = 4 * ni / np.sqrt(2)
    # k = 1.5
    # alpha = np.power(k, 4) / (1.0 + np.power(k, 4))
    # basis = gsplines.basis.Basis0101(alpha)
    # trj = gsplines.optimization.optimal_sobolev_norm(
    #     waypoints, basis, [(1, alpha), (3, 1-alpha)], execution_time)

    # trj = trj.linear_scaling_new_execution_time(2.00053)

    stop_ti = 1.29

    alpha_acc_bound = 1.98

    optimal_parametrization = opstop.minimum_time_bounded_acceleration(
        trj, stop_ti, alpha_acc_bound, str(model_file), 8)

    stop_trj = trj.compose(optimal_parametrization)

    # print(alpha)
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

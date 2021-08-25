"""
Example of how to get a curve that stops in minium time with accerelraion bouds
"""
import pathlib
import sys
import unittest
import numpy as np
import matplotlib.pyplot as plt
from .tools import debug_on

# Loads URDF model into pinocchio
from pinocchio import buildModelFromUrdf
# Stores the forward kinematics of the joints into the data argument
from pinocchio import forwardKinematics
# Updates the positions of the frames given the joints positions in data
from pinocchio import updateFramePlacements
from pinocchio import randomConfiguration
try:
    import opstop
    import pygsplines
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH))
    import opstop
    MOD_PATH_2 = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_2 = pathlib.Path(MOD_PATH_2, '..', 'build/modules/gsplines_cpp')
    sys.path.append(str(MOD_PATH_2))
    import pygsplines
    from pygsplines import optimal_sobolev_norm
    from pygsplines import BasisLegendre


def show_piecewisefunction(_q, _up_to_deriv=3, _dt=0.1, _title=''):
    dim = _q.get_codom_dim()
    fig, ax = plt.subplots(_up_to_deriv + 1, dim)
    if dim == 1:
        ax = np.array([[ax[i]] for i in range(_up_to_deriv + 1)])
    if _title:
        fig.suptitle(_title)
    t = np.arange(*_q.get_domain(), _dt)

    for i in range(0, _up_to_deriv + 1):
        q = _q.deriv(i)
        qt = q(t)
        for j in range(0, dim):
            ax[i, j].plot(t, qt[:, j])
            ax[i, j].grid()
            if i == 0:
                ax[i, j].set_title('coordinate {:d}'.format(j + 1), fontsize=8)

            if hasattr(_q, 'get_domain_breakpoints'):
                for ti in _q.get_domain_breakpoints():
                    ax[i, j].axvline(ti, alpha=0.1, color='red')

    plt.subplots_adjust(
        left=0.025,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)


def show_compare_piecewisefunction(_q, _p, _up_to_deriv=3, _dt=0.02, _title=''):
    dim = _q.get_codom_dim()
    fig, ax = plt.subplots(_up_to_deriv + 1, dim)
    if dim == 1:
        ax = np.array([[ax[i]] for i in range(_up_to_deriv + 1)])
    if _title:
        fig.suptitle(_title)
    t = np.arange(_q.get_domain()[0], _q.get_domain()[1]+_dt, _dt)
    t2 = np.arange(_p.get_domain()[0], _p.get_domain()[1]+_dt, _dt)

    for i in range(0, _up_to_deriv + 1):
        q = _q.deriv(i)
        p = _p.deriv(i)
        qt = q(t)
        pt = p(t2)
        for j in range(0, dim):
            ax[i, j].plot(t, qt[:, j], 'b')
            ax[i, j].plot(t2, pt[:, j], 'r')
            ax[i, j].grid()
            if i == 0:
                ax[i, j].set_title('coordinate {:d}'.format(j + 1), fontsize=8)

            if hasattr(_q, 'get_domain_breakpoints'):
                for ti in _q.get_domain_breakpoints():
                    ax[i, j].axvline(ti, alpha=0.1, color='red')

    plt.subplots_adjust(
        left=0.025,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)

    plt.show()


def show_compare_piecewisefunction_dk(_q, _p, _up_to_deriv=3, _dt=0.02, _title=''):

    model_file = pathlib.Path(__file__).absolute(
    ).parents[1].joinpath('urdf', 'panda_arm.urdf')

    model = buildModelFromUrdf(str(model_file))

    data = model.createData()

    dim = _q.get_codom_dim()
    fig, ax = plt.subplots(_up_to_deriv + 1, dim)
    if dim == 1:
        ax = np.array([[ax[i]] for i in range(_up_to_deriv + 1)])
    if _title:
        fig.suptitle(_title)
    t_p_spam = np.arange(_q.get_domain()[0], _q.get_domain()[1]+_dt, _dt)
    t_q_spam = np.arange(_p.get_domain()[0], _p.get_domain()[1]+_dt, _dt)

    qt = _q(t_q_spam)
    pt = _p(t_p_spam)

    joint_id = 0
    for qt_i in qt:
        forwardKinematics(model, data, qt_i)
        frame_placement = data.oMf[i]

    plt.subplots_adjust(
        left=0.025,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)

    plt.show()


class MyTest(unittest.TestCase):
    """
    """

    def __init__(self, *args, **kwargs):
        """
        """
        super(MyTest, self).__init__(*args, **kwargs)

    def error_test(self, _v_nom, _v_test):

        v_nom_inf_norm = np.max(np.abs(_v_nom))

        err_inf_norm = np.abs(_v_nom - _v_test) / v_nom_inf_norm

        assert(err_inf_norm < 1.0e-9)

    # @debug_on()
    def test(self):
        """
        """
        basis = BasisLegendre(6)
        dim = 7  # np.random.randint(1, 10)
        intervals = 4
        waypoints = np.random.rand(intervals+1, dim)
        exec_time = intervals
        trj = optimal_sobolev_norm(waypoints, basis, [(1, 3)], exec_time)

        stop_ti = 0.7*exec_time

        jerk_bound = 100.0

        diffeo = opstop.minimum_time_bouded_jerk(
            trj, stop_ti, 7*[jerk_bound])

        stop_trj = trj.compose(diffeo)

        new_exec_time = stop_trj.get_domain()[1]

        stop_trj_d1 = stop_trj.deriv()
        stop_trj_d2 = stop_trj.deriv(2)
        stop_trj_d3 = stop_trj.deriv(3)
        self.assertLess(np.max(np.abs(stop_trj_d1([new_exec_time]))), 1.0e-9)
        self.assertLess(np.max(np.abs(stop_trj_d2([new_exec_time]))), 1.0e-9)

        stop_time_span = np.arange(stop_ti, new_exec_time, 0.001)
        self.assertLess(
            np.max(np.abs(stop_trj_d3(stop_time_span))), jerk_bound*1.1)

        show_compare_piecewisefunction(trj, stop_trj)


if __name__ == '__main__':
    unittest.main()
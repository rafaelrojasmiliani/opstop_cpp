"""
Example of how to get a diffeo from the opstop
"""
import pathlib
import sys
import unittest
import numpy as np
import matplotlib.pyplot as plt
np.random.seed()
try:
    import opstop
    import gsplines
    from gsplines.optimization import optimal_sobolev_norm
    from gsplines.basis import BasisLegendre
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH))
    import opstop
    MOD_PATH_2 = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_2 = pathlib.Path(MOD_PATH_2, '..', 'build/modules/gsplines_cpp')
    sys.path.append(str(MOD_PATH_2))
    import gsplines
    from gsplines.optimization import optimal_sobolev_norm
    from gsplines.basis import BasisLegendre


def show_piecewisefunction(_q, _up_to_deriv=3, _dt=0.1, _title='', _trj_compare=None, _max_acc=None):
    dim = _q.get_codom_dim()
    fig, ax = plt.subplots(_up_to_deriv + 1, dim)
    if dim == 1:
        ax = np.array([[ax[i]] for i in range(_up_to_deriv + 1)])
    if _title:
        fig.suptitle(_title)
    t = np.arange(_q.get_domain()[0], _q.get_domain()[1]+_dt, _dt)

    for i in range(0, _up_to_deriv + 1):
        q = _q.deriv(i)
        qt = q(t)
        if _trj_compare:
            q_comp = _trj_compare.deriv(i)
            qt_comp = q_comp(t)
        for j in range(0, dim):
            ax[i, j].plot(t, qt[:, j], 'b')
            if _trj_compare:
                ax[i, j].plot(t, qt_comp[:, j], 'r')
            ax[i, j].grid()
            if i == 2:
                ax[i, j].hlines(_max_acc, t[0], t[-1])
                ax[i, j].hlines(-_max_acc, t[0], t[-1])
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


class MyTest(unittest.TestCase):
    """
    """

    def __init__(self, *args, **kwargs):
        """
        """
        super(MyTest, self).__init__(*args, **kwargs)

    def test(self):
        """
        """

        basis = BasisLegendre(6)
        dim = 7  # np.random.randint(1, 10)
        intervals = 3
        waypoints = 6*np.random.rand(intervals+1, dim) - 3
        exec_time = (intervals + 1)*0.5

        trj = optimal_sobolev_norm(waypoints, basis, [(1, 3)], exec_time)

        stop_ti = exec_time * 0.5

        xi = 0.01

        eta = 3/5*xi

        Ts = xi * (exec_time - stop_ti)
        sf = eta * (exec_time - stop_ti) + stop_ti

        diffeo = opstop.get_diffeo(stop_ti, Ts, sf)
        diffeo_diff_1 = diffeo.deriv()
        diffeo_diff_2 = diffeo_diff_1.deriv()

        time_span = np.arange(0, exec_time, 0.01)

        res = diffeo(time_span)
        res_1 = diffeo_diff_1(time_span)
        res_2 = diffeo_diff_2(time_span)

        fig, ax = plt.subplots(3, 1)
        ax[0].plot(time_span, res[:, 0])
        ax[1].plot(time_span, res_1[:, 0])
        ax[2].plot(time_span, res_2[:, 0])

        trj_stop = trj.compose(diffeo)

        a_max = np.max(np.abs(trj.deriv(2)(np.arange(0, exec_time, 0.01))))
        v_max = np.max(np.abs(trj.deriv(1)(np.arange(0, exec_time, 0.01))))

        max_acc = a_max + v_max * 16/9 / (trj.get_domain()[1] - stop_ti)/xi

        show_piecewisefunction(trj_stop, _trj_compare=trj, _max_acc=max_acc)
        plt.show()


if __name__ == '__main__':
    unittest.main()

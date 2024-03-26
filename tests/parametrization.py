"""
Test of parametrization expressions
"""
from .tools import debug_on
import pathlib
import sys
import unittest
import numpy as np
import matplotlib.pyplot as plt
import gsplines
from gsplines.functions import CanonicPolynomial
from gsplines.functions import FunctionExpression
from gsplines.functions import Identity
from gsplines.functions import ConstFunction
try:
    import opstop
except ImportError:
    import pyopstop as opstop


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
                for t_init in _q.get_domain_breakpoints():
                    ax[i, j].axvline(t_init, alpha=0.1, color='red')

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

    def error_test(self, _v_nom, _v_test):

        if np.max(np.abs(_v_nom)) < 1.0e-8:
            assert (np.max(np.abs(_v_nom-_v_test)) < 1.0e-8)

        else:
            assert (np.max(np.abs(_v_nom-_v_test)) /
                    np.max(np.abs(_v_nom)) < 1.0e-8)

    @debug_on()
    def test(self):
        """
        """
        number_of_wp = 3
        exec_time = number_of_wp - 1.0
        stop_time = exec_time * 0.9
        t_init = 0.5 * exec_time
        time_to_stop = stop_time - t_init
        s_final = t_init * (time_to_stop / t_init + 1)
        pol_coeff = [t_init, time_to_stop, 0.0,
                     -6.0 * time_to_stop + 10.0 * s_final - 10.0 * t_init,
                     8.0 * time_to_stop - 15.0 * s_final + 15.0 * t_init,
                     -3.0 * time_to_stop + 6.0 * s_final - 6.0 * t_init]
        pol = CanonicPolynomial((t_init, exec_time), pol_coeff)

        pol_2 = opstop.get_diffeo_wrt_tau(t_init, time_to_stop, s_final)

        time_span = np.arange(0, 1, 0.001)

        self.error_test(pol(time_span), pol_2(time_span))

        pol_t = pol.compose(opstop.get_tau(t_init, time_to_stop))
        pol_2_t = pol_2.compose(opstop.get_tau(t_init, time_to_stop))

        time_span = np.arange(t_init, t_init+time_to_stop, 0.001)

        self.error_test(pol_t(time_span), pol_2_t(time_span))

        diffeo = Identity((0, t_init)).concat(pol_t)
        diffeo_2 = Identity((0, t_init)).concat(pol_2_t)

        time_span = np.arange(0.0, t_init+time_to_stop, 0.1)
        self.error_test(diffeo(time_span), diffeo_2(time_span))

        for deg in range(1, 5):
            self.error_test(diffeo.deriv(deg)(time_span),
                            diffeo_2.deriv(deg)(time_span))


if __name__ == '__main__':
    unittest.main()

import pathlib
import sys
import unittest
import numpy as np
import matplotlib.pyplot as plt
try:
    from pygsplines import BasisLegendre
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH))
    import opstop
    MOD_PATH_2 = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_2 = pathlib.Path(MOD_PATH_2, '..', 'build/modules/gsplines_cpp')
    sys.path.append(str(MOD_PATH_2))
    import pygsplines


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

        stop_ti = 1.0
        exec_time = 2

        diffeo = opstop.get_diffeo(stop_ti, exec_time-stop_ti, exec_time)
        diffeo_diff_1 = diffeo.deriv()
        diffeo_diff_2 = diffeo.deriv(2)

        time_span = np.arange(0, exec_time, 0.01)

        res = diffeo(time_span)
        res_1 = diffeo_diff_1(time_span)
        res_2 = diffeo_diff_2(time_span)
        plt.plot(time_span, res[:, 0])
        plt.plot(time_span, res_1[:, 0])
        plt.show()
        diffeo_diff_1.print()


if __name__ == '__main__':
    unittest.main()

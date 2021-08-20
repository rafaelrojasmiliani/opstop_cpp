"""
Plot stuff of interest
"""


import sympy as sp
import matplotlib.pyplot as plt
import numpy as np


def main():
    xi, eta, tau = sp.symbols(r'\xi \eta \tau')
    barss = (-3*xi+6*eta-6)*tau**5 + (8*xi-15*eta + 15) * \
        tau**4 + (-6*xi+10*eta-10)*tau**3 + xi*tau + 1

    tau1 = sp.solve(barss.diff(tau, 2), tau)[2]
    tau2a = sp.solve(barss.diff(tau, 3), tau)[0]
    tau2b = sp.solve(barss.diff(tau, 3), tau)[1]

    sdot_max = sp.lambdify((xi, eta), barss.diff(
        tau).subs({tau: tau1}), 'numpy')

    sddot_a_max = sp.lambdify((xi, eta), barss.diff(
        tau, 2).subs({tau: tau2a}), 'numpy')
    sddot_b_max = sp.lambdify((xi, eta), barss.diff(
        tau, 2).subs({tau: tau2b}), 'numpy')

    dx = np.linspace(0, 4, 1000)
    dy = np.linspace(0, 4, 1000)

    xg, yg = np.meshgrid(dx, dy)

    plt.contourf(xg, yg, sdot_max(xg, yg), levels=[0, 1, 1000])

    plt.plot(dx, dx*8/15+1)

    plt.show()

    plt.contourf(xg, yg, sddot_a_max(xg, yg), levels=[0, 1, 1000])

    plt.plot(dx, dx*8/15+1)
    plt.show()

    plt.contourf(xg, yg, sddot_b_max(xg, yg), levels=[0, 1])

    plt.plot(dx, dx*8/15+1)
    plt.show()


if __name__ == "__main__":
    main()

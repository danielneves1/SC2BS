from math import *
import numpy as np
from helpfull_functions import TDMASolve

"""Units to conversion"""
tp: float = 5.391e-44  # s (1tp~5.4e-44 secs)
mp: float = 2.176e-8  #kg
lp: float = 1.616e-35
kgtoev = 5.609e35  #kg to electron volt
gigayrtosec = 1 * 60 * 60 * 24 * 365.25 * 10 ** 9 # gigayears to seconds
metertogpc = 3.24077929 * 10 ** -26 # meter to gigaparsec
c = 299_792_458  # m/s speed of light
Msun = 1.9891 * 10 ** 30  # kg mass of the sun
hbar = 1.05457182*10**-34
G = 6.6743 * 10 ** -11


def potential(x, h, psi):

    """Solves Laplacian of V = - |psi|^2 /r """

    b = np.asarray([-2 / h ** 2 for i in range(len(x))])
    b[0] = 1
    b[-1] = 1
    a = np.asarray([1 / h ** 2 for i in range(len(x))])
    a[-1] = 0
    c = np.asarray([1 /h ** 2 for i in range(len(x))])
    c[0] = 0
    norm = np.abs(psi) ** 2 / x  # u = r.psi
    #norm[0] = 0
    norm[-1] = -1
    norm[-2] = -1 + h/x[-1:]

    #  v = - self.D2inv.dot(norm)
    v = - TDMASolve(a, b, c, norm)

    return v


def der(r, dpsidr, psi, V, e, a):
       #r, x1, x2, x3, x4
    #(dx1dr =
    # dx2dr =
    # dx3dr =
    # dx4dr =)

    return (
        -(e + V/r + a/r) * psi,
          dpsidr
            )




#-------------------------------------------------------------------
calc_max = 0


def rk4Index(h, t, x1, x2, V,  e, a):
    try:
        k1 = der(t, x1, x2, V, e, a)
        k2 = der(t + h / 2, x1 + k1[0] * h / 2, x2 + k1[1] * h / 2, V, e, a)
        k3 = der(t + h / 2, x1 + k2[0] * h / 2, x2 + k2[1] * h / 2, V, e, a)
        k4 = der(t + h    , x1 + k3[0] * h,     x2 + k3[1] * h   ,  V, e, a)

        x1R = x1 + (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) * h / 6
        x2R = x2 + (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) * h / 6

        return x1R, x2R

    except:
        return inf, inf, inf, inf



def rk4(step, x_vals, uprime_initial, u_initial, v, e, alpha_ratio, **args):

    x1 = uprime_initial
    x2 = u_initial
    a = alpha_ratio

    tempData = [[x1], [x2]]


    # ver tempo
    for i in range(len(x_vals)-1):
        t = x_vals[i]
        x1, x2, = rk4Index(step, t, x1, x2, v[i], e, a)

        if calc_max != 0 and (abs(x1) > calc_max or abs(x2) > calc_max):
            break

        tempData[0].append(x1)
        tempData[1].append(x2)


    return tempData






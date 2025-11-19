
from math import *
import numpy as np
from scipy.special import genlaguerre
#from scipy.misc import derivative
from scipy.interpolate import interp1d


def hydrogen_radial_wf(r, a0, n, l):

    rho = 2*r/(a0*n)
    norm = sqrt(factorial(n-l-1) * (2/(a0*n))**3/ (2*n * factorial(n+l) ) )
    return  norm*np.exp(-rho/2) * rho**l *genlaguerre(n-l-1, 2*l+1)(rho)

def TDMASolve(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    nf = len(a)  # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))  # copy the array

    for it in range(1, nf):
        mc = ac[it] / bc[it - 1]
        bc[it] = bc[it] - mc * cc[it - 1]
        dc[it] = dc[it] - mc * dc[it - 1]

    xc = ac
    xc[-1] = dc[-1] / bc[-1]

    for il in range(nf - 2, -1, -1):
        xc[il] = (dc[il] - cc[il] * xc[il + 1]) / bc[il]

    del bc, cc, dc  # delete variables from memory

    return xc


def wf_perturbation(r, a0):
    a = 2/a0
    return r*(hydrogen_radial_wf(r, a0, 1, 0)
            + 0.23828 / a * hydrogen_radial_wf(r, a0, 2, 0)
            + 0.09891 / a * hydrogen_radial_wf(r, a0, 3, 0)
            + 0.05889 / a * hydrogen_radial_wf(r, a0, 4, 0)
            + 0.04053 / a * hydrogen_radial_wf(r, a0, 5, 0)
            + 0.03019 / a * hydrogen_radial_wf(r, a0, 6, 0)
            + 0.02366 / a * hydrogen_radial_wf(r, a0, 7, 0)
            )

#def Dperturbation(r, a0):
#    return derivative(wf_perturbation, r,1e-6, args=(a0,) )

def forward_derivative(u,i,step):
    return (-11/6 * u[i] + 3*u[i+1] -3/2 * u[i+2] +1/3 *u[i+3])/step

def xlim(slope, xn,a):
    x=xn
    n=1000
    for i in range(n):

        x1 = x - (2/a)* (1-x*a/2 -slope*np.exp(x*a/2)*np.sqrt(2/a**3))/(-2+x*a/2)

        if abs(x1-x)<1e-6:
            break

        x=x1

    if i==n-1:
        print("no convergence on xlim")

    return x


def write_to_run(x, xi, xf, xm, psi, xfinal=600, ntotal=2000, even=False, frac=10,dir="./"):

    if even:
        nlog = 0
        dx = (xfinal - xi) / (ntotal - nlog)
        xnew = np.arange(xi, xf
                              , dx)
    else:

        nlog = ntotal // frac
        dx = (xfinal - xf) / (ntotal - nlog)
        xnew = list(np.logspace(np.log(xi), np.log(xm), nlog, base=np.exp(1)))
        xnew += list(np.arange(xm + dx, xf, dx))

    interpolator = interp1d(x, psi, kind='cubic', fill_value="extrapolate")

    psinew = interpolator(xnew)



    xtoadd = list(np.arange(xf, xfinal
                            , dx))
    psinew = list(psinew) + [0.0 for i in xtoadd]
    xnew =  list(xnew) + xtoadd

    with open(dir+"WF.txt", "w") as f:
        [f.write("(" + str(i) + "," + str(0.0) + ")" + "\n") for i in psinew]

    with open(dir+"xpoints.txt", "w") as f:
        f.write(str(len(xnew)) + "\n")
        [f.write(str(i) + "\n") for i in xnew]

    return xnew, psinew

def write_input(inp, dir):
    with open(dir+"inp.txt", "w") as f:
        f.write(str(inp[0]) + "\t" + "** SpatialGridSpacing"+"\n")
        f.write(str(inp[1]) + "\t" + "** X array length" + "\n")
        f.write(str(inp[2]) + "\t" + "** xi" + "\n")
        f.write(str(inp[3]) + "\t" + "** xf" + "\n")
        f.write(str(inp[4]) + "\t" + "** fraction of points for logscale" + "\n")

        f.write("\n")

        f.write(str(inp[5]) + "\t" + "** TimeGridSpacing" + "\n")
        f.write(str(inp[6]) + "\t" + "** Max Number of time points before evap." + "\n")
        f.write(str(inp[7]) + "\t" + "** Time array length after evap" + "\n")
        f.write(str(inp[8]) + "\t" + "** initial instant" + "\n")
        f.write(str(inp[9]) + "\t" + "** final instant" + "\n")
        f.write(str(inp[10]) + "\t" + "** K" + "\n")

        f.write("\n")

        f.write(str(inp[11]) + "\t" + "** AlfaRatio" + "\n")

        f.write("\n")

        f.write(str(inp[12]) + "\t" + "** BH Initial Mass in kg" + "\n")
        f.write(str(inp[13]) + "\t" + "** BH Initial spin" + "\n")
        f.write(str(inp[14]) + "\t" + "** BH time of evaporation" + "\n")

        f.write("\n")

        f.write(str(inp[15]) + "\t" + "** Particle mass in GeV" + "\n")




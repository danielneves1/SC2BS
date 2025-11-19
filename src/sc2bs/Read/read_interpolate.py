
from scipy.interpolate import interp1d
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

"""Units to conversion"""
tp: float = 5.391e-44  # s (1tp~5.4e-44 secs)
mp: float = 2.176e-8  # kg
lp: float = 1.616e-35
kgtoev = 5.609e35  # kg to electron volt
gigayrtosec = 1 * 60 * 60 * 24 * 365.25 * 10 ** 9  # gigayears to seconds
metertogpc = 3.24077929 * 10 ** -26  # meter to gigaparsec
c = 299_792_458  # m/s speed of light
Msun = 1.9891 * 10 ** 30  # kg mass of the sun
hbar = 1.05457182 * 10 ** -34
G = 6.6743 * 10 ** -11
GeV = 1e9 / kgtoev / mp

def read_input(file_inp):
    with open(file_inp + "inp.txt") as f:
        a = f.readlines()
        intial_instant = float(a[9].split()[0])
        final_instant = float(a[10].split()[0])
        c = 3
        alfa_rat = float(a[10 + c].split()[0])
        Mbh = float(a[12 + c].split()[0]) / mp
        spin = float(a[13 + c].split()[0])
        tau = float(a[14 + c].split()[0])
        mu = float(a[16 + c].split()[0]) * GeV
    return {"alfa_rat":alfa_rat, "Mbh":Mbh, "spin":spin,"tau":tau, "mu": mu, "final_instant":final_instant,
            "intial_instant":intial_instant}

def write_input(inp):

    with open("Input/inp.txt", "w") as f:
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

def read_wf(file):
    with open(file, "r") as f:
        chunk=f.readlines()
    return [float(line.strip(" ").strip("(").strip(")\n").split(",")[0])
                + 1j*float(line.strip(" ").strip("(").strip(")\n").split(",")[1])
                for line in chunk
                if line.strip()]

def read_x(file_out):
    with open(file_out + "xpoints.txt") as f:
        a = f.readlines()

        lenx = int(len(a))
        x = np.asarray([float(a[i].split("\n")[0]) for i in range(1,lenx)])
    return x

def compute_new_alfa(t, tau, alfa_ratio_i):
    alfa_ratio = alfa_ratio_i * (1-t/tau)**(1/3)
    tau = tau * (alfa_ratio/alfa_ratio_i)**3
    return tau, alfa_ratio



def interpolate_and_write(t, tf, x, psi, xfinal, ntotal, tau,alfa_ratio_i,time_length):
    xi = x[0]
    dx = (xfinal - xi) / ntotal

    x0=int(len(x)/5)
    """This particular value DEPENDS on the CHOSEN wave function and x-grid.
     """
    xnew = np.arange(x[x0], xfinal
                     , dx)
   # plt.plot(x,np.abs(psi)**2)

    def exp_decay(x, A, B, C):
        return A * np.exp(-B * x) + C

    psireal= np.asarray(psi).real
    psiimg = np.asarray(psi).imag

    #splinereal = make_interp_spline(x, psireal, k=3)
    #freal = interp1d(x, psireal, kind='cubic', fill_value='extrapolate')
    #splineimg = make_interp_spline(x, psiimg, k=3)
    #fimg = interp1d(x, psiimg, kind='cubic', fill_value='extrapolate')

    popt, pcov = curve_fit(exp_decay, x[x0:], psireal[x0:], p0=[3, 0.5, 1])
    psinew_real = exp_decay(xnew, *popt)
    #psinew_real = freal(x)
    plt.scatter(x, psireal)
    plt.scatter(xnew, psinew_real)
    plt.show()
    popt, pcov = curve_fit(exp_decay, x[x0:], psiimg[x0:], p0=[3, 0.5, 1])
    psinew_img = exp_decay(xnew, *popt)

    plt.scatter(x, psiimg)
    plt.scatter(xnew, psinew_img)
    plt.show()


    psinew = psi[:x0]+ list(psinew_real +1j*psinew_img)

    xnew = list(x[:x0])+list(xnew)


    interpolator = interp1d(xnew, psinew, kind='linear', fill_value="extrapolate")
    xnew = np.arange(xi, xfinal
                     , dx)
    psinew = interpolator(xnew)

    tau_new, alfa_ratio = compute_new_alfa(t, tau, alfa_ratio_i)

    write_new_inps(xi=xi,xfinal=xfinal,ntotal=ntotal,frac=0,alfa_ratio=alfa_ratio, tau=tau_new,tf=tf, ti=t,
                   time_length=time_length)

    with open("Input/WF.txt", "w") as f:
        [f.write( "("+ str(i.real)+","+str(i.imag)+")" + "\n") for i in psinew]

    with open("Input/xpoints.txt", "w") as f:
        f.write(str(len(xnew)) + "\n")
        [f.write(str(i) + "\n") for i in xnew]
    return xnew, psinew

def read_t(file_out):
    with open(file_out + "tout.txt") as f:
        a = f.readlines()
        lent = int(len(a))
        t = np.asarray([float(a[i].split("\n")[0]) for i in range(lent)])
    return t



def write_new_inps(xi, xfinal, ntotal, frac, alfa_ratio,tau, tf, ti, time_length):

    inp = {"spatial eveness": "UnEven",
           "x length": ntotal,
           "xi": xi,
           "xf": xfinal,
           "fraction log pts": frac,
           "time eveness": "UnEven",
           "max time it": 5000000,
           "time array length": time_length,
           "initial time": ti,
           "final time": tf,
           "k param": 1e-4,
           "alfa_ratio": alfa_ratio,
           "BH Mass": 9999,
           "BH Spin": 9999,
           "BH evap_time": tau,
           "paricle mass in GeV": 9999,
           }
    write_input(list(inp.values()))



if __name__ == "__main__":
    """
    Here we interpolate the time evolution solution into a different x-grid.
    For our purposes we stop the simulation at alpha=1, perform this interpolation
    to then run it afterwards. This step is definitely not needed if one has enough
    computer power.
    """

    file_out="Output/"
    file_in="Input/"

    x = read_x(file_out=file_in)
    psi = read_wf(file="last_wf.txt")
    t = read_t(file_out=file_out)
    inp_dict = read_input(file_inp=file_in)

    """New values for the simulation"""
    tf = 160_000
    no_t = 60_000
    ntot = 700_000
    xf = 10_000
    x,y =interpolate_and_write(tau=inp_dict["tau"],alfa_ratio_i=inp_dict["alfa_rat"],
                          t=t[-1],tf=tf,x=x, psi=psi, xfinal=xf, ntotal=ntot,time_length=50000)



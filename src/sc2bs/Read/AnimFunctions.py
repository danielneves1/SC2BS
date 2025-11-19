from matplotlib import animation
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def read_wf(file_wf, file_x):
    with open(file_wf) as g:
        f = g.readlines()
        wf = [float(line.strip(" ").strip("(").strip(")").split(",")[0])
              + 1j * float(line.strip(" ").strip("(").strip(")\n").split(",")[1])
              for line in f
              if line.strip()]

    with open(file_x) as g:
        a = g.readlines()
        lenx = int(len(a))
        x = np.asarray([float(a[i].split("\n")[0]) for i in range(1, lenx)])

    return x, wf

################################################################################################################################
################################################################################################################################

def anim_norm(x ,tt, a, factor, s, ylim_max, results, file_out, file_in):
    def init():
        """initialize animation"""
        line2.set_data([], [])
        alf_text.set_text('')
        #  energy_text.set_text('')
        return line2, alf_text,  # energy_text

    def animate_norm(i):
        line2.set_data(x, np.abs(results[factor * i]))

        titlel.set_text(
            r'$t/\tau $=' + ' {0:1.4f}'.format(
                tt[factor * i]),)
        titler.set_text(r'$M_{BH}/M_C$' + '= {0:.2f}'.format(a[factor * i]))

        #alf_text.set_text(r'$M_{BH}/M_C$' + '= {0:.2f}'.format(a[factor * i]))
      #  alf_text2.set_text(r'$t/\tau $=' + ' {0:1.4f}'.format(tt[factor * i]),)

        return line2, alf_text

    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)


    
    ax1.set_xlim(1e-5, 1e4)
    ax1.set_xscale("log")
    titlel = ax1.set_title("", loc="left")
    titler = ax1.set_title("", loc="right")

    alf_text  = ax1.text(10 ** 2, 1.70, '')
    alf_text2 = ax1.text(10 ** 2, 1.55, '')

    line2, = ax1.plot([], [], c="blue", label="Numerical")

    cc = ["dimgrey", "grey", "silver", "pink"]
    ss=0
    # This extract the TimeIndependent solutions obtained
    for i in ["5", "2", "1"]:
        xx, yy = read_wf(file_in+"TimeIndependentSolutions/" + str(i) + "/WF.txt",
                             file_in+"TimeIndependentSolutions/" + str(i) + "/xpoints.txt")
      #  xx = np.asarray(xx) / 2
      #  yy = 2 ** 0.5 * np.asarray(yy)
        if i=="5":
            ylim_max= 1.1*max(np.abs(np.asarray(yy)))

       
        ax1.plot(xx, np.abs(np.asarray(yy)), c=cc[ss], ls="--")
       # legend.append(r"$ M_{BH} /M_{C}$="+str(i))
        ss+=1


    # this retrieves the Boson Star Solution
    xs, bs = read_wf(file_in+"TimeIndependentSolutions/BSCloud/WF.txt", file_in+"TimeIndependentSolutions/BSCloud/xpoints.txt")
   # xs = np.asarray(xs) / 2
   # bs = 2 ** 0.5 * np.asarray(bs)

    plt.plot(xs, np.abs(np.asarray(bs)), c="red", ls="-.", label="Boson Star")


    ax1.set_ylim(.0, ylim_max)

    ax1.set_ylabel(r"$\hat{r}|\hat{\psi}|$")
    ax1.set_xlabel(r"$\hat{r}$")

    anim = animation.FuncAnimation(fig, animate_norm,
                                   frames=len(tt) // factor, interval=50, init_func=init, blit=True)
    #legend += ["Boson Star"]
    ax1.legend(frameon=False)
    plt.subplots_adjust(
        left=0.15,
        right=0.90,
        top=0.90,
        bottom=0.15,
        wspace=0,
        hspace=0.0
    )
    anim.save(file_out+'Gifs/norm/u_' + str(s) + '.gif', fps=15, )
    plt.close()



def anim_norm_polar(x, tt, a, factor, s, ylim_max, results, file_out,file_in):
   
    fig, ax1 = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(6,6))
    ax1.set_title("", loc="left")
    ax1.set_title("", loc="right")
    ax1.set_yscale('log')
    ax1.set_xticks([])
    ax1.set_yticks([1],[r"$10^0$"])

   
    alf_text  = ax1.text(0, 150, '', color='white', fontsize=12)
    alf_text2 = ax1.text(0, 100, '', color='white', fontsize=12)


    n_r, n_theta = 50, 25
    r_min, r_max = 1e-4, 1e4
    r = np.logspace(np.log10(r_min), np.log10(r_max), n_r)
    theta = np.linspace(0, 2*np.pi, n_theta)
    R, Theta = np.meshgrid(r, theta)
    Z = np.zeros_like(R)

    # Initial plot
    c = ax1.pcolormesh(Theta, R, Z, shading='auto', cmap='plasma', vmin=0, vmax=1)

    # --- Initialization function ---
    def init():
        c.set_array(np.zeros(R.size))
        alf_text.set_text('')
        alf_text2.set_text('')
        return c, alf_text, alf_text2


    
    def animate_heat(i):
        r_data = np.asarray(x)
        f_r = np.asarray(np.abs(results[factor * i]))
        f = interp1d(r_data, f_r, bounds_error=False, fill_value=0)

        Z = f(R)/np.max(f(R))
        c.set_array(Z.ravel())

        #alf_text.set_text(r'$M_{BH}/M_C$ = {:.2f}'.format(a[factor * i]))
        #alf_text2.set_text(r'$t/\tau$ = {:.4f}'.format(tt[factor * i]))

        return c, alf_text, alf_text2



    anim = animation.FuncAnimation(
        fig,
        animate_heat,
        init_func=init,
        frames=max(1,len(tt)//factor),
        interval=10,
        blit=False
    )

    out_path = f"{file_out}Gifs/norm/u_polar_{s}.gif"
    anim.save(out_path, fps=15)
    plt.close()
    print(f"Saved animation to {out_path}")
################################################################################################################################
################################################################################################################################

def anim_sqrd(x, tt, a, factor, s, ylim_max, results, file_out,file_in):
    def init():
        """initialize animation"""
        line2.set_data([], [])
        alf_text.set_text('')
        #  energy_text.set_text('')
        return line2, alf_text,  # energy_text

    def animate_sqrd(i):
        line2.set_data(x, np.abs(results[factor * i]) ** 2)
        # line2.set_data(x, pot[2*i]/x)
        titlel.set_text(
            r'$t/\tau $=' + ' {0:1.4f}'.format(
                tt[factor * i]), )
        titler.set_text(r'$M_{BH}/M_C$' + '= {0:.2f}'.format(a[factor * i]))
        #   alf_text.set_text(r'$ M_{BH}/M_C $= {0:1.2f}'.format(a[factor * i]))

        return line2, alf_text

    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)


    ax1.set_xlim(1e-5, 1e4)
    ax1.set_xscale("log")
    # ax1.set_yscale("log")
    titlel = ax1.set_title("", loc="left")
    titler = ax1.set_title("", loc="right")
    alf_text = ax1.text(10, 0.4, '')
    #line1, = ax1.plot([], [], c="orange")
    line2, = ax1.plot([], [],c="blue", label= "Numerical")


    cc = ["dimgrey", "grey", "silver", "pink"]
    ss = 0
    for i in ["5", "2", "1"]:

        xx, yy = read_wf(file_in+"TimeIndependentSolutions/" + str(i) + "/WF.txt",
                             file_in+"TimeIndependentSolutions/" + str(i) + "/xpoints.txt")
      #  xx = np.asarray(xx) / 2
      #  yy = 2 ** 0.5 * np.asarray(yy)
        if i == "5":
            ylim_max = 1.1 * max(np.abs(np.asarray(yy))**2)

        ax1.plot(xx, np.abs(np.asarray(yy))**2, c=cc[ss], ls="--")

        ss += 1
    #this retrieves the Boson Star Solution
    xs, bs = read_wf(file_in+"TimeIndependentSolutions/BSCloud/WF.txt", file_in+"TimeIndependentSolutions/BSCloud/xpoints.txt")
  #  xs = np.asarray(xs) / 2
  #  bs = 2 ** 0.5 * np.asarray(bs)
    plt.plot(xs, np.abs(np.asarray(bs))**2, c="red", ls="-." , label="Boson Star")

    ax1.legend(frameon=False)

    ax1.set_ylabel(r"$\hat{r}|\hat{\psi}|$")
    ax1.set_xlabel(r"$\hat{r}$")
    ax1.set_ylim(.0, ylim_max)

    anim = animation.FuncAnimation(fig, animate_sqrd,
                                   frames=len(tt) // factor, interval=10, init_func=init, blit=True)
    plt.subplots_adjust(
        left=0.15,
        right=0.90,
        top=0.90,
        bottom=0.15,
        wspace=0,
        hspace=0.0
    )
    anim.save(file_out+'Gifs/sqrd/u2_' + str(s) + '.gif', fps=15, )
    
    plt.close()

#######################################################################################################################
################################################# ZOOM ################################################################
#######################################################################################################################

def anim_sqrd_Zoom(x, tt, a,  factor, dd, ylim_max, results, file_out, file_in):
    def init():
        """initialize animation"""
        line2.set_data([], [])

        #  energy_text.set_text('')
        return line2,
    def animate_sqrd(i):
        line2.set_data(x, np.abs(results[factor * i]) ** 2)
        # line2.set_data(x, pot[2*i]/x)
        #title.set_text('t = {0:1.2f}'.format(tt[factor * i]))
        titlel.set_text(
            r'$t/\tau $=' + ' {0:1.4f}'.format(
                tt[factor * i]), )
        titler.set_text(r'$M_{BH}/M_C$' + '= {0:.2f}'.format(a[factor * i]))
        return line2,

    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)
    ax1.set_xlim(1, x[-1])
    titlel = ax1.set_title("", loc="left")
    titler = ax1.set_title("", loc="right")
    line2, = ax1.plot([], [], c="blue", label="Numerical")

    cc = ["dimgrey", "grey", "silver", "pink"]

    ss=0
    # This extract the TimeIndependent solutions obtained
    for i in ["5", "2", "1"]:
        xx, yy = read_wf(file_in+"TimeIndependentSolutions/" + str(i) + "/WF.txt",
                         file_in+"TimeIndependentSolutions/" + str(i) + "/xpoints.txt")
      #  xx = np.asarray(xx) / 2
      #  yy = 2 ** 0.5 * np.asarray(yy)
        if i == "5":
            ylim_max = 1.1 * max(np.abs(np.asarray(yy))**2)
        ax1.plot(xx, np.abs(np.asarray(yy)) ** 2, c=cc[ss], ls="--")

        ss += 1
    # this retrieves the Boson Star Solution
    xs, bs = read_wf(file_in+"TimeIndependentSolutions/BSCloud/WF.txt",file_in+ "TimeIndependentSolutions/BSCloud/xpoints.txt")
    #xs = np.asarray(xs) / 2
   # bs = 2 ** 0.5 * np.asarray(bs)
    plt.plot(xs, np.abs(np.asarray(bs)) ** 2, c="red", ls="-.", label="Boson Star")

    ax1.legend(frameon=False)
    ax1.set_ylim(.0,ylim_max)
    ax1.set_xlim(.0, 25)
    ax1.set_ylabel(r"$\hat{r}|\hat{\psi}|$")
    ax1.set_xlabel(r"$\hat{r}$")
    anim = animation.FuncAnimation(fig, animate_sqrd,
                                   frames=int((len(tt)) / factor), interval=100, init_func=init, blit=True)

    anim.save(file_out+'Gifs/sqrd/u2_Zoom_' + str(dd) + '.gif', fps=15, )
    plt.close()

################################################################################################################################
################################################################################################################################

def anim_norm_Zoom(x, tt, a, factor, dd, ylim_max, results, file_out, file_in ):
    def init():
        """initialize animation"""
        line2.set_data([], [])

        #  energy_text.set_text('')
        return line2,
    def animate_norm(i):
        line2.set_data(x, np.abs(results[factor * i]) )
        # line2.set_data(x, pot[2*i]/x)
       # title.set_text('t = {0:1.2f}'.format(tt[factor * i]))
        titlel.set_text(
            r'$t/\tau $=' + ' {0:1.4f}'.format(
                tt[factor * i]), )
        titler.set_text(r'$M_{BH}/M_C$' + '= {0:.2f}'.format(a[factor * i]))
        return line2,

    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)
    ax1.set_xlim(1, x[-1])
    titlel = ax1.set_title("", loc="left")
    titler = ax1.set_title("", loc="right")
    line2, = ax1.plot([], [], c="blue", label="Numerical")

    cc = ["dimgrey", "grey", "silver", "pink"]

    ss = 0
    #This extract the TimeIndependent solutions obtained
    for i in ["5", "2", "1"]:
        xx, yy = read_wf(file_in+"TimeIndependentSolutions/" + str(i) + "/WF.txt",
                         file_in+"TimeIndependentSolutions/" + str(i) + "/xpoints.txt")
        #xx = np.asarray(xx) / 2
       # yy = 2 ** 0.5 * np.asarray(yy)
        if i == "5":
            ylim_max = 1.1 * max(np.abs(np.asarray(yy)))
        ax1.plot(xx, np.abs(np.asarray(yy)) , c=cc[ss], ls="--")

        ss += 1

    #this retrieves the Boson Star Solution
    xs, bs = read_wf(file_in+"TimeIndependentSolutions/BSCloud/WF.txt", file_in+"TimeIndependentSolutions/BSCloud/xpoints.txt")
    #xs = np.asarray(xs) / 2
    #bs = 2 ** 0.5 * np.asarray(bs)
   
    plt.plot(xs, np.abs(np.asarray(bs)) , c="red", ls="-.", label ="Boson Star")

    ax1.legend(frameon=False)

    ax1.set_ylim(.0, ylim_max)
    ax1.set_xlim(.0, 25)
    ax1.set_ylabel(r"$\hat{r}|\hat{\psi}|$")
    ax1.set_xlabel(r"$\hat{r}$")
    anim = animation.FuncAnimation(fig, animate_norm,
                                   frames=int((len(tt)) / factor), interval=100, init_func=init, blit=True)

    anim.save(file_out+'Gifs/norm/u_Zoom_' + str(dd) + '.gif', fps=15, )
    plt.close()





import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from AnimFunctions import read_wf


def side_by_side_plot(x, y_data, t_data,n_cols,n_rows, m_data,file_in,file_out):
    font = {
        "family": "Times New Roman",
        'size': 40}

    matplotlib.rc('font', **font)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 30 ), sharex=True, sharey=True, )
    plt.subplots_adjust(left=0.16, right=0.95, bottom=0.2, top=0.99)

    axes = axes.flatten()
    cc = ["dimgrey", "grey", "silver", "pink"]

    for i, ax in enumerate(axes):


        ax.tick_params(
            axis='both',
            which='major',
            length=10,
            width=4,
            direction='inout',
            color='black',
            labelsize=40
        )


        ss = 0
        for j in ["5", "2", "1"]:
            xx, yy = read_wf(file_in+"TimeIndependentSolutions/" + str(j) + "/WF.txt",
                                        file_in+ "TimeIndependentSolutions/" + str(j) + "/xpoints.txt")
            #xx = np.asarray(xx) / 2
            #yy = 2 ** 0.5 * np.asarray(yy)
            if j == "5":
                ylim_max = 1.1 * max(np.abs(np.asarray(yy)))

            ax.plot(xx, np.abs(np.asarray(yy)), c=cc[ss], ls="--", linewidth=2.5)

            ss += 1
        xs, bs = read_wf(file_in+ "TimeIndependentSolutions/BSCloud/WF.txt", file_in+ "TimeIndependentSolutions/BSCloud/xpoints.txt")

        ax.plot(xs, np.abs(np.asarray(bs)), c="red", ls="--", label="Boson Star", linewidth=2.5)
        
        ax.set_ylabel(r"$\hat{r}|\hat{\psi}|$")
        ax.set_xlabel(r"$\hat{r}$")

        ax.plot(x[i], y_data[i] , c="blue", linewidth=4)
        
        alf_text =  ax.text(0.15 * 10 ** 1, 1.10, '')
        alf_text2 = ax.text(0.15 * 10 ** 1, 0.95, '')
        alf_text.set_text(r'$M_{BH}/M_c $' + '= {0:.2f}'.format(m_data[i]))
        alf_text2.set_text(r'$t/\tau $=' + ' {0:1.4f}'.format(t_data[i]), )
        ax.set_xlim(1e-3, 2e2)
        ax.set_ylim(0,  ylim_max)

        ax.set_xscale("log")
        # This hides  y-axis labels except for the leftmost plots
        if i % 2 != 0:  # not in first column
            ax.set_ylabel('')
            ax.tick_params(labelleft=False)

        # This hides the x-axis labels except for bottom row
        if i < 5 :  # top row
            ax.set_xlabel('')
            ax.tick_params(labelbottom=False)
        if i == 7:
            ax.set_xticks([1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2])  # positions
            ax.set_xticklabels(["", r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$'])

    plt.subplots_adjust(
        left=0.1,
        right=0.98,
        top=0.99,
        bottom=0.05,
        wspace=0,
        hspace=0.0
    )

    plt.savefig(file_out+"SidebySideFrames.pdf")


from scipy.interpolate import interp1d

def side_by_side_polar_plot(r_datas, f_rs):
    font = {
        'weight': 'bold',
        'size': 30}

    matplotlib.rc('font', **font)
    fig, axes = plt.subplots(4, 2, figsize=(20, 30),
                             subplot_kw={'projection': 'polar'}, sharex=True, sharey=True)

    axes = axes.flatten()
    cc = ["dimgrey", "grey", "silver", "pink"]

    for i, ax_polar in enumerate(axes):

        r_data = np.asarray(r_datas[i])
        f_r = np.asarray(f_rs[i])
        f_r = f_r
        f = interp1d(r_data, f_r, bounds_error=False, fill_value=0)

        # Polar grid (logarithmic radius)
        n_r = 200
        n_theta = 200
        r_min, r_max = 1e-2, 200

        r = np.logspace(np.log10(r_min), np.log10(r_max), n_r)
        theta = np.linspace(0, 2 * np.pi, n_theta)

        R, Theta = np.meshgrid(r, theta)
        Z = f(R)


        c = ax_polar.pcolormesh(Theta, R, Z, shading='auto', cmap='plasma', vmax=1.6, vmin=0)

        ax_polar.set_yscale('log')
        for label in ax_polar.get_yticklabels():
            label.set_color('white')
            label.set_fontsize(20)

        ax_polar.set_xticks([])  # remove angular ticks
        ax_polar.set_yticks([])

    plt.subplots_adjust(
        left=0.0,
        right=0.98,
        top=0.95,
        bottom=0.08,
        wspace=-0.5,
        hspace=0.0
    )
    cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.02])
    fig.colorbar(c, cax=cbar_ax, orientation='horizontal')


if __name__=="__main__":

    #if we already saved the specific results we just need to load them and plot afterward
    x1 = np.load("x1.npy")
    x2 = np.load("x2.npy")

    res1 = np.load("res1.npy")
    res2 = np.load("res2.npy")

    mass1 = np.load("mass1.npy")
    mass2 = np.load("mass2.npy")

    time1 = np.load("time1.npy")
    time2 = np.load("time2.npy")

    x = [list(x1) for i in range(len(res1))] + [list(x2) for i in range(len(res2))]
    res = list(res1)+ list(res2)
    time = list(time1)+ list(time2)
    mass = list(mass1)+ list(mass2)

    side_by_side_polar_plot(x, res)
    side_by_side_plot(x=x,y_data=res, t_data=time, m_data=mass)


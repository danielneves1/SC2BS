import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def plot_energy(x, f, file_out):
    font = {
        "family": "Times New Roman",
        'size': 15}

    matplotlib.rc('font', **font)

    fig, ax = plt.subplots(figsize=(6, 3))
    fig.tight_layout()
    fig.subplots_adjust(left=0.16, right=0.95, top=0.99)
    ax.plot(x, f, c="blue", label="Numerical")

   
    ax.set_xlim(x[0], x[-1])

  #  ax.set_yticks([-0.5, -0.4, -0.3, -0.2, -0.1, 0])
    
    ax.set_ylabel(r"$\langle \hat{E} \rangle $")

    x = np.linspace(-1, max(x), 100)
    ax.plot(x, [0 for i in x], c="black")
    ax.plot(x, [-0.027 for i in x], c="red", ls="--", label="Boson Star ")

    ax.set_xlabel(r"$t / \tau $")

    plt.savefig(file_out + "Energy.pdf")

    plt.close()
    print("Energy Plot Done")


def plot_norm(x, f, file_out):
    plt.plot(x, 1-np.asarray(f))
    plt.ylabel(r"$norm shift$")

    plt.title("1 - Norm")
    plt.xlabel(r"$time$")

    plt.savefig(file_out+"Norm.pdf")
    plt.close()


if __name__ == "__main__":

    """ Simply plot the energy and norm"""

    filedir = "Bounded simulation/"
    file_out = filedir + "Output/"
    file_in =   filedir +"Input/"

    from readers import read_instants, read_norm, read_energy,read_input

    print("Reading Inputs")
    inp_dict = read_input(file_in)


    ttt = read_instants(file_out)
    t0 = ttt[0]
    ttt = np.asarray(ttt / (inp_dict["tau"] + t0))  # normalize in terms of evaporation time


    norm = read_norm(file_out)
    energy = read_energy(file_out)

    plot_energy(ttt, energy, file_out)
    plot_norm(ttt, norm, file_out)


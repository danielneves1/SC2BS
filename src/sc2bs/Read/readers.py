import numpy as np
from tqdm import tqdm

def read_instants(file_out):
    with open(file_out + "instants.txt") as f:
        a = f.readlines()

        lent = int(len(a))
        t = np.asarray([float(a[i].split("\n")[0]) for i in range(lent)])
    return t

def read_input(file_inp):
    with open(file_inp + "inp.txt") as f:
        a = f.readlines()
        intial_instant = float(a[9].split()[0])
        final_instant = float(a[10].split()[0])
        c = 3
        alfa_rat = float(a[10 + c].split()[0])
        Mbh = float(a[12 + c].split()[0])
        spin = float(a[13 + c].split()[0])
        tau = float(a[14 + c].split()[0])
        mu = float(a[16 + c].split()[0])
    return {"alfa_rat":alfa_rat, "Mbh":Mbh, "spin":spin,"tau":tau, "mu": mu, "final_instant":final_instant,
            "intial_instant":intial_instant}

def read_norm(file_out):
    with open(file_out + "norm.txt") as f:
        a = f.readlines()

        lenx = int(len(a))
        x = np.asarray([float(a[i].split("\n")[0]) for i in range(lenx)])
    return x



def read_energy(file_out):
    with open(file_out + "energy.txt") as f:
        a = f.readlines()
        lenx = int(len(a))
        x = np.asarray([float(a[i].split("\n")[0]) for i in range(lenx)])
    return x

def read_x(file_out, reduced):
    if reduced:
        st = "xpoints_reduced.txt"
    else:
        st = "xout.txt"
    with open(file_out + st) as f:
        a = f.readlines()

        lenx = int(len(a))
        x = np.asarray([float(a[i].split("\n")[0]) for i in range(lenx)])
    return x

def read_t(file_out):
    with open(file_out + "tout.txt") as f:
        a = f.readlines()

        lent = int(len(a))
        t = np.asarray([float(a[i].split("\n")[0]) for i in range(lent)])
    return t



# Function to convert chunk data to floats
def parse_chunk(chunk, f, x):

    return [float(line.strip(" ").strip("(").strip(")").split(",")[0])
                + 1j*float(line.strip(" ").strip("(").strip(")").split(",")[1])
                for line in tqdm(chunk.splitlines())
                if line.strip()]


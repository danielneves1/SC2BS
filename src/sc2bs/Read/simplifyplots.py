import concurrent.futures

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')
from tqdm import tqdm
from readers import read_x, read_t, read_input, parse_chunk


import os

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

# Function to read a chunk of the file

val=0
s=0
dummy2=0
dd = 0
flag = False

def read_and_reduce(futures, x, extra, missing, total, file_out, numb,all_indices):

    lenx=len(x)

    global s
    global dummy2
    global flag


    for f in futures:

        f = f.result()
        if f!=[]:

            to_add = extra + list(f[:missing])
            if to_add == []:
                s=1
            else:
                wf = to_add
                wf = np.asarray(wf)
                psinew = wf[all_indices]

                # write in the same format as in the f90 files/ Crank-Nicolson.exe output

                with open(file_out + "WF_reduced.txt", "a") as sss:
                    [sss.write("(" + str(np.real(i)) + "," + str(np.imag(i)) + ")" + "\n") for i in psinew]

            dummy = (len(f) - missing) // lenx
            remainder = (len(f) - missing) % lenx

            print("Shortening to" + str(numb) + " pts" + str(s + 1) + "/" + str(total))
            for j in range(dummy):
                psinew = np.asarray([f[missing + i + j * lenx] for i in range(lenx)])[all_indices]
                with open(file_out + "WF_reduced.txt", "a") as sss:

                    [sss.write("(" + str(np.real(i)) + "," + str(np.imag(i)) + ")" + "\n") for i in psinew]

            extra = list(f[missing + dummy * lenx:])
            missing = lenx - remainder


            dummy2 += dummy


            f.clear()
            s+=1

    return extra, missing

def read_chunk(file, start, size, file_size, file_name):
    global val
    file.seek(val)

    chunk = file.read(size)
    val = file.tell()
    if not chunk.endswith('\n') and start + size < file_size:

     chunk += file.readline()
     val = file.tell()
    return chunk

# Main function to process the file
def read_file_in_parallel_to_reduce(file_in, file_out, file_name, num_workers, runs):
    file_name1 = file_name
    file_name = file_out + file_name
    file_size = os.path.getsize(file_name)
    #chunk_size = file_size // num_workers //runs

    print("Reading Inputs")
    inp_dict = read_input(file_in)
    x = read_x(file_out, False)
    t = read_t(file_out)

    print("Reading Wave Functions")
    extra = []
    missing = 0
    start = 0
    finnish = file_size // runs

    global val
    val = start

    #initial and final points of the new x-grid

    ini = int(np.log10(x[0]))
    fin = int(np.log10(x[-1]))

    numb = 200 #Total number of points to reduce to.

    all_indices = []

    for i in range(ini, fin + 1):

        cond = (x > 10 ** i) & (x < 10 ** (i + 1))
        valid_indices = np.where(cond)[0]

        if len(valid_indices) >= numb:

            selected = np.linspace(0, len(valid_indices) - 1, num=numb, dtype=int)
            sampled_indices = valid_indices[selected]

            all_indices += list(sampled_indices)
        else:
            all_indices += list(valid_indices)

    xnew = x[all_indices] #new x-grid

    #write new x-grid
    with open(file_out+"xpoints_reduced.txt", "w") as sss:

        [sss.write(str(i) + "\n") for i in xnew]

    chunk_size = (file_size - start )//num_workers//runs

    for i in range(runs):

        with open(file_name, 'r') as file:

            chunks = [read_chunk(file, ini, chunk_size, file_size, file_name1) for ini in range(start, finnish, chunk_size)]

        # Process chunks in parallel
        # Process chunks in parallel using ProcessPoolExecutor
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(parse_chunk, chunk, file_name1, x) for chunk in tqdm(chunks)]



            extra, missing = read_and_reduce(futures, x, extra, missing, int(num_workers*runs), file_out, numb, all_indices)
        futures = []
        start = finnish
        finnish += finnish

    return dd
# Usage example:
if __name__ == "__main__":
    import argparse
   
    parser = argparse.ArgumentParser()
    parser.add_argument("-workers", "--workers", help="Number of Cores")
    parser.add_argument("-runs", "--runs", help="Number of times the file is divided")
    parser.add_argument("-input", "--input", help="Input file location")
    parser.add_argument("-output", "--output", help="Output file location")

    #manual input:

    # filedir = "Bounded simulation/"
    # file_out = filedir + "Output/"
    # file_in =   filedir +"Input/"
    # runs = 1
    # workers = 1

    args = parser.parse_args()
    workers=int(args.workers)
    runs = int(args.runs)

    file_out = args.output
    file_in = args.input

    print("Simplifying / Reducing data")
    print("")
    if os.path.exists(file_out + "WF_reduced.txt"):
    
        print("`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´`´ ")
        print("Another WF_reduced.txt exists at ", f"{file_out}")
        print("This code will write above said file")
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print(" ")
        print("###################################")
        print("############ ENDING ###############")
        print("###################################")
        raise Exception
        
    read_file_in_parallel_to_reduce(file_in=file_in,file_out=file_out, file_name="wfout.txt", num_workers=workers, runs=runs)



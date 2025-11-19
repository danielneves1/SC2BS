import concurrent.futures

import numpy as np
import matplotlib

matplotlib.use('Agg')
from tqdm import tqdm

from AnimFunctions import (anim_norm, anim_sqrd, anim_sqrd_Zoom,
                           anim_norm_Zoom, anim_norm_polar)

from SideSide import side_by_side_plot

from readers import (read_instants, read_x, read_t, read_norm,
                     read_energy, read_input, parse_chunk)


from BasicPlots import plot_energy, plot_norm

from simplifyplots import read_file_in_parallel_to_reduce

import os


val = 0
s=0
dummy2=0
dd = 0
flag = False

def WF_gif(futures, x, t, inp_dict, extra, missing, total, file_out, file_in):
    """Here we process the wave function data and produce Gifs"""
    #Notice that x and the wf (futures) are written in Ruffini Bonazzola units, thus different
    #normalizations may be needed depending on what we wish to animate


    lenx=len(x)

    #x=x/2 #different units
    x = x
    factor = 1
    tau= inp_dict["tau"]
    alfa_rat= inp_dict["alfa_rat"]
    t0 = inp_dict["intial_instant"]
    tf = inp_dict["final_instant"]
    taup = t0 + tau

    global s
    global dummy2
    global flag

    for f in futures:

        f = f.result()
        if f!=[]:
            case=0
           # to_add = [extra + list(2**0.5*np.asarray(f[:missing]))] #these are for different unit system
            to_add = [extra + list(np.asarray(f[:missing]))]

            dummy =  (len(f) - missing) // lenx - case
            remainder = (len(f) - missing) % lenx

            # we only wish to save each wave function such that
            # result[i] corresponds to the wave function at the instant t[i].
            # It may occur, when reading the file in parallel,since we are dividing the file
            # into multiple chunks, that we do not obtain the full wave function (for example, we stop at the middle ).
            # To circunvent this issue, we store this partial wave function array (array called extra)
            # and add it to the next data chunk.

            result = [ np.asarray([ f[missing + i + j * lenx]
                              for i in range(lenx)])
                for j in range(case,dummy)]

            # result = [2 ** 0.5 * np.asarray([f[missing + i + j * lenx]
             #                          for i in range(lenx)])
              #  for j in range(case,dummy)]




            tt = [t[dummy2 + i] for i in range(case, dummy)]

            extra = list(np.asarray(f[missing + dummy * lenx:]))
            #extra = list(2**0.5*np.asarray(f[missing + dummy * lenx:]))


            missing = lenx - remainder

            if to_add==[[]]:
                results = np.asarray(result)
                a = np.asarray([alfa_rat * (1 - (i - t0) / tau) ** (1 / 3) if i - t0 <= tau else 0 for i in tt])


                timetoplot = np.asarray([ i/taup for i in tt])
            else:
                results = np.asarray(to_add + result)
                tt.append(t[dummy2 + dummy ])
                a = np.asarray([alfa_rat * (1 - (i - t0) / tau) ** (1 / 3) if i - t0 <= tau else 0 for i in tt])
                timetoplot = np.asarray([i / taup for i in tt])

            dummy2 += dummy
            global ylim_max, ylim_max_z
            global ylim2_max, ylim2_max_z
            if s==0:

               # ylim_max = 1.1 * max(abs(results[0])) # If you don't know the exact value of the initial Wave function
               # ylim2_max = 1.1 * max(abs(results[0])**2)
               #here I am fixing the y_lim to a particular value when I compute the animations(inside the fucntion)
               ylim_max = 99999
               ylim2_max = 99999


            print("Producing Animations " + str(s + 1) + "/" + str(total))

            # 2D polar-plot of |\psi|
            print("Gif: WaveFunction Abs (Polar)")
            anim_norm_polar(x, timetoplot, a, factor, s, ylim_max, results, file_out, file_in)

           # 1D plot of |\psi|
            print("Gif: WaveFunction Abs")
            anim_norm(x, timetoplot, a, factor, s, ylim_max, results, file_out, file_in)

            # 1D plot of |\psi|^2
            print("Gif: WaveFunction Squared")
            anim_sqrd(x, timetoplot, a, factor, s, ylim_max, results, file_out, file_in)
           
            

            #If we wish to save specific wf's
            # such that no longer we need to read the entire file, we can simply save it
            # using:
            if len(futures)==1:
                n_cols=2
                n_rows=4
                #We will plot as follows for the 4x2
                #might need to change frame indices
                #indxs =[ 0, 10,
                #        20, 30,
                #        40, 50,
                #        60, 80]
                
                indxs= [ int(k/(n_cols*n_rows) * len(results)) for k in range(int(n_cols*n_rows))] # example
                
                #a1 = np.asarray([np.abs(results[i]) for i in indxs])
                #np.save("res2.npy", a1)
                #a1 =  np.asarray([timetoplot[i] for i in indxs])
                #np.save("time2.npy", a1)
                #a1 = np.asarray([a[i] for i in indxs])
                #np.save("mass2.npy", a1)
                #np.save("x2.npy", x)

                print("Done")
                #Plots side by side 4x2 images from at the instants labeled in indxs
                
                side_by_side_plot(
                x = [x for i in indxs],
                y_data= [np.abs(results[i]) for i in indxs],
                t_data= np.asarray([timetoplot[i] for i in indxs]),
                m_data=np.asarray([a[i] for i in indxs]),
                file_in=file_in, file_out=file_out, n_cols=n_cols, n_rows=n_rows)




            global dd

            tg=[i for i in range(len(tt)) if tau<=tt[i]]
            if tg != []:
                idk=1
                flag=True

            if flag:

                factor=1
                if idk == 1:
                    start= tt.index(tt[tg[0]])
                    idk=2
                else:
                    start=0

                tt = tt[start:]
                results = results[start:]
                if dd == 0:
                   # ylim_max_z = 1.1 * max(abs(results[0]))
                   # ylim2_max_z = 1.1 * max(abs(results[0]) ** 2)
                    ylim_max_z = 99999
                    ylim2_max_z = 99999
                    ini=tt[0]
                    
                #In case one wishes to plot linearly and not in log-scale uncomment the following lines
               # anim_norm_Zoom(x, timetoplot, a,  factor, dd,  ylim_max_z, results, file_out,file_in)
               # anim_sqrd_Zoom(x, timetoplot, a,  factor, dd, ylim2_max_z, results, file_out, file_in)
          
                dd+=1
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

   # if file_name=="wfout.txt":

    #    if not chunk.endswith('\n') and start + size < file_size:

     #       chunk += file.readline()
      #      val = file.tell()

   # elif file_name == "Potentialout.txt":

    #    if not chunk.endswith('\n') and start + size < file_size:
    #        chunk += file.readline()
    #        val = file.tell()
            # Complete the line in case it's split
    return chunk


# Main function to process the file
def read_file_in_parallel(file_in, file_out, file_name, num_workers, runs,reduced):
    file_name1= file_name
    file_name = file_out + file_name
    file_size = os.path.getsize(file_name)


    print("Reading Inputs")
    inp_dict = read_input(file_in)
    x = read_x(file_out, reduced)

    # These are the time steps for each outputed wave-function.
    t = read_t(file_out)

    #These are the total time steps used in the simulation
    ttt = read_instants(file_out)
    t0 = ttt[0]
    ttt = np.asarray(ttt / (inp_dict["tau"] + t0)) #normalize in terms of evaporation time

    norm = read_norm(file_out)
    energy = read_energy(file_out)

    plot_energy(ttt, energy, file_out)
    plot_norm(ttt, norm, file_out)

    print("Reading Wave Functions")

    extra = []
    missing = 0
    start = 0
    finnish = file_size // runs

    global val
    val = start

    chunk_size = (file_size - start )//num_workers//runs

    for i in range(runs):

        with open(file_name, 'r') as file:

            chunks = [read_chunk(file, ini, chunk_size, file_size, file_name1) for ini in range(start, finnish, chunk_size)]

        # Process chunks in parallel
        # Process chunks in parallel using ProcessPoolExecutor
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(parse_chunk, chunk, file_name1, x) for chunk in tqdm(chunks)]



            extra, missing = WF_gif(futures, x, t, inp_dict, extra, missing, int(num_workers*runs), file_out, file_in)
        futures = [] #can delete from memory if needed
        start = finnish
        finnish += finnish

    return dd

# Usage example:
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-workers", "--workers", help="Number of Cores")
    parser.add_argument("-runs", "--runs", help="Number of times the file is divided into N Cores")
    parser.add_argument("-input", "--input", help="Input file location")
    parser.add_argument("-output", "--output", help="Output file location")
    parser.add_argument("-reduced", "--reduced", help="Reduced output or not")

    #This code allows reading the data in parallel using a specific
    # number of cores, named workers, and further allows to divide the reading
    # into different number of runs, named runs.
    # This essentially divides the file into workers * runs pieces.

    # "manual input":

    #filedir = "Bounded simulation/"
    #reduced = True
    #file_out = filedir + "Output/"
    #file_in =   filedir +"Input/"
    #runs = 1
    #workers = 1

    args = parser.parse_args()
    workers=int(args.workers)
    runs = int(args.runs)

    file_out = args.output
    file_in = args.input
    reduced = args.reduced
    try:
        os.mkdir(file_out + "Gifs/norm")
    except FileExistsError:

        pass
        
    try:
        os.mkdir(file_out + "Gifs/sqrd")
    except FileExistsError:

        pass
        
    if reduced == "True" or reduced == "true":
        reduced= True
    elif reduced == "False" or reduced == "false":
        reduced= False
    else:
        print(f"{reduced=}")
        print("And it should be either True or False")
        raise Exception
   # reduced = True
   
    
    #This mainly concerns to the fact that,
    # if there are too many points in the simulation
    # it may be helpful to previously reduce the number of points
    # via an interpolation or particularly choosing a few, for example.
    # This will expedite the production of the animation/plots without loosing much information.
    # For visual purposes we believe it suffices. Here we called it WF_reduced.txt .

    #If you wish to reduce the number of points to plot
    # you can use this function:

    #read_file_in_parallel_to_reduce(file_in=file_in,file_out=file_out, file_name="wfout.txt", num_workers=workers, runs=runs)
    # before the next piece of code, using reduced = True.

    # Also please consider using
    #
    # num_workers = 1
    # runs = 1
    #
    # since now the reduced data may be small enough to read it in one go,
    # and there will be no need to join the multiple Gifs afterwards


    if reduced:

        n_zooms=read_file_in_parallel(file_in=file_in,file_out=file_out, file_name="WF_reduced.txt", num_workers=workers,
                                      runs=runs, reduced=reduced)

    else:
        n_zooms = read_file_in_parallel(file_in=file_in, file_out=file_out, file_name="wfout.txt",
                                        num_workers=workers,
                                        runs=runs, reduced=reduced)


    #If ran in parallel, this code will give multiple gifs which can be concatenated using the following:

    #For the case where the Wave Function is reduced and sufficient memory is available we can skip this step
    # and simply use one worker and one run.

    #from JoinGifs import concatenate_gifs

    #n_gifs = int(runs * workers)
    #giflocation = "/Output/Gifs/norm/"
    #gifname = "u_" #in this example each gif is labelled as u_0.gif, u_1.gif, u_2.gif etc...
    #gif_paths = [giflocation + gifname + str(s)+".gif" for s in range(0, n_gifs)]# Paths to the gifs you wish to join

    #output_path = giflocation + "u.gif"  # Path to save the combined gif
    #concatenate_gifs(gif_paths, output_path)


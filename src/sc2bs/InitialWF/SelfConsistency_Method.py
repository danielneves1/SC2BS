import matplotlib.pyplot as plt
import math
from RK4 import rk4, potential
import numpy as np
from helpfull_functions import (wf_perturbation,
                                forward_derivative, write_to_run,
                                xlim, write_input)


def energy_cycle(energy, xm, args):
    energy_flag= False
    wf_flag = False
    no_convergence_flag = False
    if energy > 0:
        raise Exception
   # l = list(str(energy))
    # ["-", #1, #2, ..., ".", #d1, #d2, ... ]
    # start energy variation on the same order
   # if abs(energy) >= 1:
   #     order_magnitude = l.index(".") - 1

    #    de = 10 ** (order_magnitude - 1)
   # else:
   #     ini = l.index(".")
   #     for i in range(ini + 1, len(l)):
    #        if l[i] == "0":
    #            continue
     #       else:
     #           de = 10 ** - (i - 2)
     
    xmax=max(x_vals)
    n = 5
    dn = 0.5
    order_magnitude = math.floor(math.log10(abs(energy)))
    de = 10 ** (order_magnitude )
    tol_abs=1e-5
    tol_rel = 1e-8
    #print(de)
    for j in range(20):
        

        """These cycles attempt to find the correct energy value
        such that at the end of the curve does not go to either 
        positive or negative infinity. Thus, looking for the energies 
        where the curve rises to +infinity and where it decreases towards 
        -infinity, we get that the correct energy value must lie between 
        these two solutions. Thus, obtaining the correct solution.
         """
  
        for i in range(12):
           # print(energy, de)
            dudr, u = rk4(e=energy,**args)
            
            maxu = max(u)
            val = [ dudr[k] for k in range(len(dudr)) if  abs(dudr[k]) < tol_abs + tol_rel*abs(maxu)  and x_vals[k] > xm]
           # print(min(np.abs(dudr)), tol_abs + tol_rel*abs(maxu))
            val_idx = [k for k in range(len(dudr)) if abs(dudr[k]) < tol_abs + tol_rel*abs(maxu) and x_vals[k] > xm]
            du_s = np.asarray(dudr)[(np.asarray(x_vals) > xm) & (np.asarray(x_vals) < xmax )]
            no_zeros = not(np.any((du_s[:-1] < 0) & (du_s[1:] > 0)) )
           # plt.plot(x_vals,u, label=f"{i=}")
            if u[-1]>0 and i==0:
                print("ERR: Possibly energy is too high")
                print("")
                print("Modifiying initial energy guess")
                print("...")
            if u[-1] > 0:
                break
            
            #check if it has zeros (maybe energy is too high and computes higher modes)
           # signs = np.sign([u[i] for i in range(len(u)) if x_vals[i] < 0.7 * xmax ])
        #    no_zeros = not(np.any(signs[:-1] * signs[1:] < 0))
        
            energy_old= energy
            energy -= de
            #no_zeros=True
            
            if abs(energy_old-energy)<1e-14 and val!=[] and no_zeros:
                energy_flag=True
                #no_convergence_flag=False
                break
            if n - dn <=1:
                n=1
            else:
                n-=dn
          #  print(energy)
       # plt.ylim(-2,10)
       # plt.legend()
       # plt.show()

        if energy_flag:
            print("ENERGY CONVERGED: Check curve manually (n might be to low)")

            min_deriv = min(val, key=abs)
            ind = val_idx[val.index(min_deriv)]
            print("Min of WF is ", min(u), "Min of derivative", min_deriv)
            u = np.asarray(u)
            #val=[i for i in range(len(dudr)) if -0.001 > dudr[i] <0.01 and x_vals[i]>xm]

            u[ind:]=0


            break
        energy += de
        de = 0.1 * de
        #print(no_zeros)
       
        #if val !=[]:
        
          #  print("Might check Inflection point")
    if j == 19:
        no_convergence_flag = True
        #print("Couldn't Converge Wave Function")
   # plt.plot(x_vals,u)
    #plt.show()
    if not(no_zeros) and not(energy_flag):
        print("Solution can have zeros.")
    return energy, u , no_convergence_flag , dudr


def convergence(old, new):
    """Basic convergence example test"""
    total=len(old)
    avg_diff = sum([abs(old[i]-new[i]) for i in range(total)])/total

    if avg_diff<1e-6:
        return True
    else:
        return False




def guess(r):

    return wf_perturbation(r,a0=2/alpha_ratio)

def drivativ(r):
    return Dperturbation(r,a0=2/alpha_ratio)


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

if __name__=="__main__":

    """ Here we compute the Initial Conditions/Initial Wave-Function (WF) suited to our
      problem in a self-consistent manner. We also, define the spatial grid 
      for the CrankNicolson Simulation. Furthermore, we write the Input file 
      and set the inputs to later evolve it using the Crank-Nicolson method.
      """

    import argparse



    parser = argparse.ArgumentParser()
    parser.add_argument("-alpha_ratio", "--alpha_ratio", help="alpha_ratio")
    parser.add_argument("-tau", "--tau", help="Evaporation Time of BH")
    parser.add_argument("-tinitial", "--tinitial", help="Initial time of simulation")
    parser.add_argument("-tfinal", "--tfinal", help="Final time of simulation")
    parser.add_argument("-time_steps", "--time_steps", help="time steps after BH evaporation ")
    
    parser.add_argument("-ntotal", "--ntotal", help= "Number of x-grid points")
    parser.add_argument("-xfinal", "--xfinal", help= "Final point in x-grid")
    parser.add_argument("-xinitial", "--xinitial", help="Initial point in x-grid")
    parser.add_argument("-dir", "--dir", help="Directory to write Input")
    
    alpha_ratio = float(parser.parse_args().alpha_ratio)
    ntotal = int(parser.parse_args().ntotal)
    xfinal = float(parser.parse_args().xfinal)
    xinitial = float(parser.parse_args().xinitial)
    dir=parser.parse_args().dir
    time_steps= int(parser.parse_args().time_steps) #notice that these time steps are only set after t=tau thus dt=(tfinal-tau)/time_steps
    tinitial= float(parser.parse_args().tinitial)
    tfinal= float(parser.parse_args().tfinal)
    tau=float(parser.parse_args().tau)
    
    
    """
    First Simply define a x-grid to converge the initial solution. 
    Only after we interpolate/extrapolate it until the inputed xfinal.
    For alpha ratios of order O(1) to O(10) xf = 30 is alright. 
    If alpha_ratio> 10 or 20 consider change both the step and a smaller xf.
    """
    if alpha_ratio < 1:
        xf = 70
    if 10>alpha_ratio >= 1:
        xf = 30
    if alpha_ratio>=10:
        xf = 10
        
    step = 0.01
    xi = xinitial
    x_vals = np.arange(xi, xf + step
                       , step)

    frac = 99999999 # legacy value


    old_u = guess(x_vals)
    u_initial = old_u[0]
    uprime_initial = forward_derivative(old_u, 0,step)


    xinflection = 3 * 4 / alpha_ratio

    energy_0 = - 1/4*alpha_ratio**2 #- 5/16*alpha_ratio


    no = 1
    found_n = False

    """ Convergence Cycle """
    for k in range(1, 100):

        v = potential(x=x_vals,h=step,psi=old_u)
        plt.plot(x_vals,old_u)
        rk4args = {"step": step, "x_vals": x_vals, "uprime_initial":uprime_initial,
                   "u_initial":u_initial,"v":v , "alpha_ratio": alpha_ratio }#this exact order
        no_conv_flag = True

        if k==1:
            for no in range(1,6):
                if no_conv_flag:
                    energy = energy_0*no
                    energy, u, no_conv_flag, dudr = energy_cycle(energy, xinflection, args= rk4args)
                    if no_conv_flag:
                        print("Trying different starting energy, attempt:", no)
                    else:
                        nn=no
        else:

            energy = energy_0 * nn
            energy, u, no_conv_flag, dudr = energy_cycle(energy, xinflection, args=rk4args)
        if no==6 or no_conv_flag==True:
            print("!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!!!FAILED TO CONVERGE!!!")
            print("!!!!!!!!!!!!!!!!!!!!!!!!")
            raise Exception

        print("Iteration:", k," energy:", energy, "Norm=", sum(u**2*step))
        u = u / np.sqrt(sum(u ** 2 * step))
        u_initial=u[0]
        uprime_initial=forward_derivative(u,0,step)

        """Plot to see the evolution of the solution"""
     #   plt.plot(x_vals, u)

        if convergence(old_u, u):
            print("########################################")
            print("##### SOLUTION CONVERGED at it:", k,"####")
            print("########################################")
            break
        else:
            old_u = u

  #  plt.show()

    """Verify solution manually through plots"""
   # v = potential(x=x_vals,h=step,psi=u)
   # plt.scatter(x_vals, v)
   # plt.show()

    #plt.scatter(x_vals, u)
    #plt.show()
    plt.close()
    xinflection = 3 * 4/alpha_ratio
    xm = xlim(-0.01,xinflection, alpha_ratio)


    x,y = write_to_run(x=x_vals,xi=xi,xf=xf,xm=xm,
                 psi=u, xfinal=xfinal,ntotal=ntotal,even=True,frac=frac, dir=dir)

    inp = {"spatial eveness" :"UnEven",
           "x length": ntotal,
           "xi":xi,
           "xf":xfinal,
           "fraction log pts": frac,
            "time eveness" :"UnEven",
           "max time it": 50000000,
           "time array length": time_steps,
           "initial time":  tinitial,
           "final time": tfinal,
           "k param": 1e-4,
           "alpha_ratio": alpha_ratio,
           "BH Mass": 9999,
           "BH Spin": 9999,
           "BH evap_time": tau,
           "paricle mass in GeV": 9999,
    }
    write_input(list(inp.values()), dir)

    """Check if Number of Points accurately describes the Solution at least visually."""
    
    print("Check if the exponential part has enough points.")

    plt.title("Initial State Scatter Plot")
    plt.scatter(x, y)
    plt.xscale("log")
    plt.xlim(xi/10, 4*xinflection)
    plt.savefig("InitialWF.png")



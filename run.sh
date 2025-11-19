#!/bin/bash


cd src/sc2bs/
# Set number of cores
file_out="./CrankNicolson/Output/"
file_in="./CrankNicolson/Input/"


#x-grid properties
xinitial=1e-3
xfinal=100
ntotal=4000

#black hole + field initial conditions
alpha_ratio=0.5  # Black hole to cloud mass ratio
tau=20 #evaporation time

#time grid properties
tinitial=0
tfinal=100 #tfinal should be greater than tau
time_steps=400 #notice that these time steps are only set after t=tau thus dt=(tfinal-tau)/time_steps
workers=1
runs=1

reduced=False
#Here we set reduced to False because the code that reads and plots
# by default assumes no reduction/simplification, thus
# it would output error if no reduced file exists.


#generate initial WaveFunction to run and input parameters for the CrankNicolson
echo " "
echo "Generating Initial Wave Function using input values in run.sh file"
echo " "
python3 InitialWF/SelfConsistency_Method.py -alpha_ratio="$alpha_ratio" -ntotal="$ntotal" -xinitial="$xinitial" -xfinal="$xfinal" -dir="$file_in" -tau="$tau" -tinitial="$tinitial" -tfinal="$tfinal" -time_steps="$time_steps"

echo " "
echo "Check InitialWaveFunction.png file."
echo " "
#go to CrankNicolson Folder to generate executable
cd CrankNicolson

echo " "
echo "Evolving Initial Wave Function Through Crank-Nicolson"
echo " "
gfortran -o CrankNicolson.exe variables.f90 cloud_potential.f90 ini.f90 Main.f90 math_methods.f90 predictor.f90 wave_function.f90 writers.f90 instants.f90  modulation.f90 -ffree-line-length-none

#Run executable
./CrankNicolson.exe

#read data, do plots and Gifs.
cd ..
############################################################################################
 #If we wish to reduce/simplify the data we can choose only a few points to plot
 #If not please comment the following lines
echo " "
echo " "
echo "Simplifying Data for smaller computer power" "----->"  "$reduced"
echo " "

if [ "$reduced" = True ]; then

python3 Read/simplifyplots.py -workers="$workers" -runs="$runs" -input="$file_in" -output="$file_out"
echo " "
echo "Reading and Plotting Data"
echo " "
#read data, do plots and Gifs.
python3 Read/Main_Reader.py -workers="$workers" -runs="$runs" -input="$file_in" -output="$file_out" -reduced="$reduced"

else
echo " "
echo "Reading and Plotting Data"
echo " "
python3 Read/Main_Reader.py -workers="$workers" -runs="$runs" -input="$file_in" -output="$file_out" -reduced="$reduced"
#The file is called "WF_reduced.txt" and its associated x-grid "xpoints_reduced.txt"
#This particular name is required if reduction since it is not an input.(can easily be changed)
############################################################################################
fi

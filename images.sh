#!/bin/bash

#script to simply obtain the images
# Set number of cores
cd src/sc2bs/
file_out="./CrankNicolson/Output/"
file_in="./CrankNicolson/Input/"

reduced=False

workers=1
runs=1


############################################################################################
 #If we wish to reduce/simplify the data we can choose only a few points to plot
 #If not please comment the following lines
#python3 Read/simplifyplots.py -workers="$workers" -runs="$runs" -input="$file_in" -output="$file_out"
#reduced=True

#The file is called "WF_reduced.txt" and its associated x-grid "xpoints_reduced.txt"
#This particular name is required if reduction since it is not an input.(can easily be changed)
############################################################################################

#read data, do plots and Gifs.
python3 Read/Main_Reader.py -workers="$workers" -runs="$runs" -input="$file_in" -output="$file_out" -reduced="$reduced"

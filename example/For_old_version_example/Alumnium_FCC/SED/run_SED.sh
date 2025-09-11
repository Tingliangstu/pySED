#!/bin/bash
# ********************************************************************************************
# ******* One need to edit the input_SED.in file and make this parameters correctly **********
# ********************************************************************************************

# Firstly, we need to run calculate SED mode and get the SED data when make sure the parameters in input_SED.in file is 
# totolly set

plot_SED=0               # Set to 0 to calculate SED, run SED mode
sed  -i  "41s/1/${plot_SED}/g" input_SED.in 

#conda activate pySED     # Activate the conda inv (maybe one dont need it)

pySED input_SED.in

# Secondly, plot SED and fitting mode, one maybe need to run this mode multiple times

plot_SED=1
sed  -i  "41s/0/${plot_SED}/g" input_SED.in

pySED input_SED.in


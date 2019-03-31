#!/bin/bash

awk '{print $1}' E_tot.txt > time_aux
awk '{print $2}' E_kinpot.txt > pot_aux
awk '{print $3}' E_kinpot.txt > kin_aux
awk '{print $2}' E_tot.txt > tot_aux
awk '{print $2}' PressTemp.txt > press_aux
awk '{print $3}' PressTemp.txt > temp_aux 
awk '{print $2}' Momentum_z.txt > moment_aux
paste Momentum_xy.txt moment_aux > Momentum.txt
paste  time_aux pot_aux kin_aux tot_aux press_aux temp_aux > Results.txt
rm *_aux

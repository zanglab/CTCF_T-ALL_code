#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BTN3A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BTN3A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTN3A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTNL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTNL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTNL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BTNL8 

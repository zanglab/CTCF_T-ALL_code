#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EXOSC4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOSC9 

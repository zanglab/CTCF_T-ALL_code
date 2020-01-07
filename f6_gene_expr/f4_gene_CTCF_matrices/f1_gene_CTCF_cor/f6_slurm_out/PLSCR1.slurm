#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLSCR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLSCR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLSCR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLSCR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLSCR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLSCR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLTP 

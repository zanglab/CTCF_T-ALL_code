#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NME6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NME6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NME7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NME8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NME9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMNAT1 

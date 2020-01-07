#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NMNAT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NMNAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMNAT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMRAL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMRK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMRK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NMT1 

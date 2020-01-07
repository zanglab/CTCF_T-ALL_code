#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NR1I3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NR1I3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR2C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR2C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR2C2AP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR2E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR2E3 

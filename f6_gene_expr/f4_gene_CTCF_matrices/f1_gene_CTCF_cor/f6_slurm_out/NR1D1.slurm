#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NR1D1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NR1D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR1D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR1H2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR1H3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR1H4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NR1I2 

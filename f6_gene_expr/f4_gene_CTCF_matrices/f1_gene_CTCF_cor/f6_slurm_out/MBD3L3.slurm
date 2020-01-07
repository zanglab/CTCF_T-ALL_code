#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MBD3L3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MBD3L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBL1P 

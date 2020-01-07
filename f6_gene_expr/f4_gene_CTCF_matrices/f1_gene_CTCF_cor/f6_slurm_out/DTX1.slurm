#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DTX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DTX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTX3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTYMK 

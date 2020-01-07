#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MPPED2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MPPED2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPRIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPV17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPV17L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPV17L2 

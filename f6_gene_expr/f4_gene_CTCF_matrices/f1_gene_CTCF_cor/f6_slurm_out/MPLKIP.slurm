#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MPLKIP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MPLKIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPND 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MPP3 

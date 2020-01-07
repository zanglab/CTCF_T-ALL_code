#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLFN5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLFN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLFNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLIRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLIT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLIT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLIT3 

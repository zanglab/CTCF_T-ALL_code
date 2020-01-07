#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAR2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAR2P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FARP1 

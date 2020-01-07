#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPLP0 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPLP0P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPLP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPLP2 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPGRIP1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPGRIP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPH3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPH3AL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPIA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL10A 

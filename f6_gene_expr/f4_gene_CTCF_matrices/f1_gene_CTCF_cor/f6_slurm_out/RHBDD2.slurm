#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RHBDD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RHBDD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHBDD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHBDF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHBDF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHBDL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RHBDL2 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PFKM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PFKM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFKP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFN1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PFN3 

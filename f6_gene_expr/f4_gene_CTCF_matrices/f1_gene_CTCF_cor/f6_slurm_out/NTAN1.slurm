#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NTAN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NTAN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTHL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NTMT1 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTPRZ1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTPRZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTRH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTRH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTRHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTTG1 

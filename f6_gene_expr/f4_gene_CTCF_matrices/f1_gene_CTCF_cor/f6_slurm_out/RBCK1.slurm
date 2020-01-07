#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBCK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBCK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBFA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBFOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBFOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBFOX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBKS 

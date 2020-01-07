#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OAZ1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OAZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OAZ2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OAZ3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OBP2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OBP2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OBSCN 

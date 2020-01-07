#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIBF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIBF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PICALM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PICK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PID1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIDD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIEZO1 

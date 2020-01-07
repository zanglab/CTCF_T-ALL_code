#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMED5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED7-TICAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED9 

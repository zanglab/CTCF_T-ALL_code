#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLAUR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLAUR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLBD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLBD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCB2 

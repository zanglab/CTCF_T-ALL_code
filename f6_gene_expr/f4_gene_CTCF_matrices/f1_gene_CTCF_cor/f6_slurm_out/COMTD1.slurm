#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o COMTD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i COMTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COPG1 

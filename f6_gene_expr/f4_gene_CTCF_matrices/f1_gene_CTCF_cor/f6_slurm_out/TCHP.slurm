#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TCHP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TCHP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCIRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCL1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCL1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCN2 

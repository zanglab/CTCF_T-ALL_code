#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPPL2C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPPL2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPPL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRED1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRED2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPRED3 

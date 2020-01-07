#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OPALIN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OPALIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPCML 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPHN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPLAH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPN1LW 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OPN1MW 

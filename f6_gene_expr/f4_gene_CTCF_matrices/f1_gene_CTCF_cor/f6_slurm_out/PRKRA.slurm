#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRKRA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKRIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKXP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRKY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRL 

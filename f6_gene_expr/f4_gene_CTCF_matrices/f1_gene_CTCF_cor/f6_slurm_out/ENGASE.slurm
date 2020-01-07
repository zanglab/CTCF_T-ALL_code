#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ENGASE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ENGASE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENHO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENKD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENKUR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENO2 

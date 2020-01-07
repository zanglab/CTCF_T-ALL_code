#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DUXA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DUXA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUXAP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DVL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DVL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DVL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DXO 

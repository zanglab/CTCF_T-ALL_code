#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP5S.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP5S 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6AP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6AP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6AP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0A2 

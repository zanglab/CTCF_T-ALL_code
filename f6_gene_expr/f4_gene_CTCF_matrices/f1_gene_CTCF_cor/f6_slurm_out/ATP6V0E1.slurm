#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP6V0E1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0E1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0E2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V1B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V1B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V1C1 

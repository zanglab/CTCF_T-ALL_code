#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP6V0A4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0CP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V0D2 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP6V1G2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V1G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V1G2-DDX39B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V1G3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP6V1H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP7B 

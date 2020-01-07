#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CACNA1D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA1F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA1G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA1H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNA1I 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o JADE2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i JADE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JADE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JAG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JAG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JAGN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JAK1 

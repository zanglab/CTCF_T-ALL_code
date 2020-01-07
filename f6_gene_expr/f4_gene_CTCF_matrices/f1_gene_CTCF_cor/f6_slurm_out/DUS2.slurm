#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DUS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DUS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUS3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUS4L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DUSP11 

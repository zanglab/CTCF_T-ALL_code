#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDI2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDIAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDIT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDIT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDIT4L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDN 

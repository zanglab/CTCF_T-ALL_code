#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTRNR2L2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTRNR2L7 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTBP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTCH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTCP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTDH 

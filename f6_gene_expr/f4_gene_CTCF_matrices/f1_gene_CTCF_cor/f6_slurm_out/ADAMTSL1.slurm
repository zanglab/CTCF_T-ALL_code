#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADAMTSL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTSL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTSL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTSL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTSL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAP1 

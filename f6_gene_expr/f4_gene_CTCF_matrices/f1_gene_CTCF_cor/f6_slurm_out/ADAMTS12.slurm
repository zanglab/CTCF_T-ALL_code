#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADAMTS12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAMTS17 

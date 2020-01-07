#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNNI1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNNI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNNI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNNI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNNI3K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNNT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNNT2 

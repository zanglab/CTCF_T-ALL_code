#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCNO.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNYL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNYL2 

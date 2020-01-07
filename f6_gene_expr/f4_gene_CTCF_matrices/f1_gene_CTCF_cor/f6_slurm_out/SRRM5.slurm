#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRRM5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRRM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRRT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRSF12 

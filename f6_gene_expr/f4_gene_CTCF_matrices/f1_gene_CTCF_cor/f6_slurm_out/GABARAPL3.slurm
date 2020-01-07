#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GABARAPL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GABARAPL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABBR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABBR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABPB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABPB2 

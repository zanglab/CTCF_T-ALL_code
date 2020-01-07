#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CACNG4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNG6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACNG8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CACTIN 

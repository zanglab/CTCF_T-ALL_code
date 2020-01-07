#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHA5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHAC1 

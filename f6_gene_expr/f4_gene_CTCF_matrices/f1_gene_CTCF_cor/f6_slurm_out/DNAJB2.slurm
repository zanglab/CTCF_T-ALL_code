#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAJB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJB7 

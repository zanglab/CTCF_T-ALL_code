#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAJB8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC12 

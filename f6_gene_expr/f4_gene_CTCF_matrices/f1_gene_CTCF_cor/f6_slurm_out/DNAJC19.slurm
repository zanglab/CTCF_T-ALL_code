#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAJC19.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC25 

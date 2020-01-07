#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAJC25-GNG10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC25-GNG10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC4 

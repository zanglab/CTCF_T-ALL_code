#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAJC13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC18 

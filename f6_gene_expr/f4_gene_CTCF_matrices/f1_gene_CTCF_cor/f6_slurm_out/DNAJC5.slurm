#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNAJC5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC5G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNAJC8 

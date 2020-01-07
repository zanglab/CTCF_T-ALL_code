#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GTDC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GTDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2A1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GTF2E1 

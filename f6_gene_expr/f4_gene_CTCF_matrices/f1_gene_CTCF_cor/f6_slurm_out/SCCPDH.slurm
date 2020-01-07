#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCCPDH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCCPDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCEL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCFD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCFD2 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SHOX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SHOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHPK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHPRH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHQ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHROOM1 

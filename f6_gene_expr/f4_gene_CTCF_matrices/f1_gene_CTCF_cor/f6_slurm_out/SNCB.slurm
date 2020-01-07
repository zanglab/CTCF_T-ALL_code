#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNCB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNCB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNCG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNED1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNIP1 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BBS5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BBS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBS9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BCAN 

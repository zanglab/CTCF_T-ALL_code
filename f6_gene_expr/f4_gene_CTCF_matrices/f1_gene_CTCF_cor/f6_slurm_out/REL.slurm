#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o REL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i REL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RELA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RELB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RELL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RELL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RELN 

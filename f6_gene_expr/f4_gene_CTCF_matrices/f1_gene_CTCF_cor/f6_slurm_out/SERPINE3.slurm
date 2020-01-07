#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SERPINE3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPING1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINI1 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1GALT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1GALT1C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1GALT1C1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1QB 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LAMP5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMTOR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMTOR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMTOR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMTOR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMTOR5 

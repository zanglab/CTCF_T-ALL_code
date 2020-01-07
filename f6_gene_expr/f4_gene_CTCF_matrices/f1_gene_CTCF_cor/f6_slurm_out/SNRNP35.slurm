#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNRNP35.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRNP35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRNP40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRNP48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRNP70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPA1 

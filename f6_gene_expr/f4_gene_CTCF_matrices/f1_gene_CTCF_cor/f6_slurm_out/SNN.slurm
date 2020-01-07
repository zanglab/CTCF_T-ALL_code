#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRNP200 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRNP25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRNP27 

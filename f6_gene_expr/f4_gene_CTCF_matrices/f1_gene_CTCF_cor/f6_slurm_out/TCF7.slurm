#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TCF7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF7L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF7L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCFL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCHH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCHHL1 

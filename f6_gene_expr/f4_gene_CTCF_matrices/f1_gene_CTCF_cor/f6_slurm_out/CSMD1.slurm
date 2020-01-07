#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CSMD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CSMD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSMD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSMD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSN1S1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSN3 

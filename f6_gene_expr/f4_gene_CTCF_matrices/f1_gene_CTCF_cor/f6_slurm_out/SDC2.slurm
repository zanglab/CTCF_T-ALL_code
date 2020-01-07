#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SDC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDC4P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDCBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SDCBP2 

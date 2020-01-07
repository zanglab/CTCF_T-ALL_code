#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CSGALNACT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CSGALNACT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSGALNACT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSHL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSK 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDUFA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA4L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA7 

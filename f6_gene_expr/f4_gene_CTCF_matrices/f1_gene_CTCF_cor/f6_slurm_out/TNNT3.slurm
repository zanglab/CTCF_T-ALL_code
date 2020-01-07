#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNNT3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNNT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNPO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNPO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNPO3 

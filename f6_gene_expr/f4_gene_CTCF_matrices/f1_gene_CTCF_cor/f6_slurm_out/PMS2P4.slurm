#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PMS2P4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS2P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS2P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS2P7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS2P9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMVK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNCK 

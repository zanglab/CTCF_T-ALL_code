#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PMS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS2CL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS2P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMS2P3 

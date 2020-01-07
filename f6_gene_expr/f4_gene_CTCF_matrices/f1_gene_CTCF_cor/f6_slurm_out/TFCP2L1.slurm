#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TFCP2L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TFCP2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFDP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFDP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFDP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TFEB 

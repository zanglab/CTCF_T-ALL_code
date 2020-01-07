#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UGDH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UGDH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGGT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGGT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UGT1A10 

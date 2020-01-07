#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DCP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DCP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCSTAMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCT 

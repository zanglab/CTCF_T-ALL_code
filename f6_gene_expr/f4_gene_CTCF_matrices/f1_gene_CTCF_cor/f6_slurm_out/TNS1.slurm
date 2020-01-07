#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNXA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNXB 

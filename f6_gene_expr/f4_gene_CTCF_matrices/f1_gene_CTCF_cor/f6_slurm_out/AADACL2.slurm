#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AADACL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AADACL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AADACL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AADACL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AADACP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AADAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AAED1 

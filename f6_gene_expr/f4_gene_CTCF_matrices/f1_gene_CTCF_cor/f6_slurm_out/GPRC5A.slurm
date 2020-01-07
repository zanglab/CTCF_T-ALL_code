#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPRC5A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRC5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRC5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRC5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRC5D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRC6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRIN1 

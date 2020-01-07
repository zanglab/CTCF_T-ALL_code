#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPIRE2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIRE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPNS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPNS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPNS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPO11 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBR4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBR5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBR7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBTD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBTF 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o YTHDF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i YTHDF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YWHAB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YWHAE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YWHAEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YWHAEP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i YWHAG 

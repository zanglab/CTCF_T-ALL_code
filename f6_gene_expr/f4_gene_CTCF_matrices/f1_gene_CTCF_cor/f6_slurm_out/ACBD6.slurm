#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACBD6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACBD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACBD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACCS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACCSL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACE 

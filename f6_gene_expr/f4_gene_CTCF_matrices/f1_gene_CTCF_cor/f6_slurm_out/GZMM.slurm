#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GZMM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GZMM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H1F0 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H1FNT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H1FOO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H1FX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i H2AFB1 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APEX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APEX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APH1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APH1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i API5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APLF 

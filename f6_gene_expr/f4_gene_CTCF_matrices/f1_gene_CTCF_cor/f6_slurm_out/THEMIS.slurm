#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o THEMIS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i THEMIS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THEMIS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THG1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THNSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THNSL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i THOC1 

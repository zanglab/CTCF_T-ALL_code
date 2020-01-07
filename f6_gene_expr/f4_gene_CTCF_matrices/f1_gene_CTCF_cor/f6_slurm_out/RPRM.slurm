#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPRM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPRM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPRML 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS10-NUDT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS10P7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS11 

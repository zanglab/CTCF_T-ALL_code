#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RELT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RELT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RENBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i REP15 

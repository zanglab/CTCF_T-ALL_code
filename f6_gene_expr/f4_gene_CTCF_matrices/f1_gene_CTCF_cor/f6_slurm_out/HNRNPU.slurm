#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HNRNPU.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPUL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPUL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HNRNPUL2-BSCL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOGA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOMER1 

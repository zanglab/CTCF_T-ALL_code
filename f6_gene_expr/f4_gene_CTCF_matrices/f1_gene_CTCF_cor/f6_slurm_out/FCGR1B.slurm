#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FCGR1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FCGR1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCGR2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCGR2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCGR2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCGR3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FCGR3B 

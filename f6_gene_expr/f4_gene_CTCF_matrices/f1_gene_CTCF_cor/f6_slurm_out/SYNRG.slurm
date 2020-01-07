#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYNRG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNRG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYS1-DBNDD2 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CNTNAP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTNAP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTNAP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNTROB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i COA3 

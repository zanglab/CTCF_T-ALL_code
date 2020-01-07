#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o QPCT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i QPCT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QPCTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QPRT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QRFP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QRFPR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i QRICH1 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IPP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IPP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPPK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCA1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IQCC 

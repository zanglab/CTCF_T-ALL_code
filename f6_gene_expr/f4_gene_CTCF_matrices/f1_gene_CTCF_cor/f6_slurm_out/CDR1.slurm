#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDR2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDRT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDRT15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDRT15L2 

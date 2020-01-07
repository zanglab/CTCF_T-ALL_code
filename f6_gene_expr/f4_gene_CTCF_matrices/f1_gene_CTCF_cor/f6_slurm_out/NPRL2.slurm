#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPRL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPRL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPRL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPSR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPTX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPTX2 

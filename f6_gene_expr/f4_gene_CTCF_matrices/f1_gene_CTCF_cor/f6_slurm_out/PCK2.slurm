#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCLO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCMTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCMTD2 

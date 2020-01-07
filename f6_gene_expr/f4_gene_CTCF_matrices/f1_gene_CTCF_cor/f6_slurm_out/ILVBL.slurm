#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ILVBL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ILVBL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMMP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMMP2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IMP4 

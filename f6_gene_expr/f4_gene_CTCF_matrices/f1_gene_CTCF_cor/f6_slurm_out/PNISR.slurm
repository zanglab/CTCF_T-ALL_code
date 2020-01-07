#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PNISR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PNISR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNKD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNKP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNLDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNLIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNLIPRP1 

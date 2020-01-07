#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UNK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UNK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNKL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UPF2 

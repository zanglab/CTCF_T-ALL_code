#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GOSR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GOSR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOT1L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GOT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GP1BA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GP1BB 

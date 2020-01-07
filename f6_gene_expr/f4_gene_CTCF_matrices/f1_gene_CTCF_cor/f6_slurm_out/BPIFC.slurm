#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BPIFC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPNT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPTF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRAF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRAP 

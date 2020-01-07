#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAV2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAV2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAV3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBARP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBFA2T2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBFA2T3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBFB 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RTCA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RTCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTCB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTEL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTEL1-TNFRSF6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RTKN 

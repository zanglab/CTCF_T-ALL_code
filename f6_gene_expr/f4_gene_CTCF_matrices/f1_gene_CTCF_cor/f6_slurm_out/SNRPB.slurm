#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNRPB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPD2P2 

#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNRPD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNRPN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNTA1 

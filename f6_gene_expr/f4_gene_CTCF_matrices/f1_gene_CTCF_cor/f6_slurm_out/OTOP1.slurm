#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OTOP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTOS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OTP 

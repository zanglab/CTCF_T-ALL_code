#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IL17F.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IL17F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL17RA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL17RB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL17RC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL17RD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IL17RE 

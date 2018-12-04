#!/bin/bash 
# 
#$ -cwd 
#$ -V 
#$ -j y 
#$ -S /bin/bash 
# 
 
R --vanilla <./code/mainScript.R> mainScript.Rout

#!/bin/sh
export EMAIL=prlarsen@anl.gov

qsub -n 64 -t 02:00:00 -A HLRedshift -q prod -M ${EMAIL} ./subme.sh


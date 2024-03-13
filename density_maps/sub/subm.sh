#!/bin/sh
export EMAIL=prlarsen@anl.gov

qsub -I -n 1 -q debug -t 01:00:00 -A LastJourney -M ${EMAIL}
qsub -n 32 -t 12:00:00 -A LastJourney -q default -M ${EMAIL} ./Cobalt_cooley.sh


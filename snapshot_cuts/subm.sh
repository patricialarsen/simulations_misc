#!/bin/sh
export EMAIL=prlarsen@anl.gov

qsub -n 4 -t 0:30:00 -A LastJourney -q debug -M ${EMAIL} ./Cobalt_cooley.sh
#qsub  -I -n 1 -q debug -t 00:40:00 -A LastJourney -M ${EMAIL} 
#./Cobalt_cooley.sh

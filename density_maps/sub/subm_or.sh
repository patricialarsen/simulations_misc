#!/bin/sh
export EMAIL=prlarsen@anl.gov

#qsub -I -n 10 -q debug -t 02:00:00 -A LastJourney -M ${EMAIL}
#qsub -n 16 -t 1:00:00 -A LastJourney -q debug -M ${EMAIL} ./Cobalt_cooley_or.sh
#qsub -n 32 -t 6:00:00 -A LastJourney -q default -M ${EMAIL} ./Cobalt_cooley_or_32.sh
#qsub -n 4 -t 1:00:00 -A LastJourney -q debug -M ${EMAIL} ./Cobalt_cooley_or_downs.sh

qsub -n 8 -t 0:40:00 -A LastJourney -q debug -M ${EMAIL} ./Cobalt_cooley_bc.sh


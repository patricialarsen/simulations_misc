#!/bin/sh

export EMAIL=prlarsen@anl.gov
#qsub -q default -n 32 -t 120 -A LastJourney -M $EMAIL ./Cobalt_cooley.sh

#qsub -q default -n 1 -t 240 -A LastJourney -M $EMAIL ./Cobalt_cooley.sh
qsub -q default -n 1 -t 30 -A LastJourney -M $EMAIL ./Cobalt_cooley.sh

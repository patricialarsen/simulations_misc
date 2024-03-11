#!/bin/sh

export EMAIL=prlarsen@anl.gov
#qsub -q debug -n 4 -t 120 -A LastJourney -I
#qsub -q debug -n 8 -t 60 -A LastJourney -M $EMAIL -I 
#qsub -q debug -n 8 -t 20 -A LastJourney -M $EMAIL ./cobalt_pl.sh

qsub -q debug -n 4 -t 120 -A LastJourney -M $EMAIL ./cobalt_pl_2023.sh


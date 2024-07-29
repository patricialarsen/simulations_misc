#!/bin/sh
export EMAIL=prlarsen@anl.gov

qsub -n 32 -t 12:00:00 -A LastJourney -q default -M ${EMAIL} ./Cobalt_polaris.sh


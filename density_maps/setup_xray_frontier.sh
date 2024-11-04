#!/bin/bash
source /ccs/home/prlarsen/codes/HACC/env/bashrc.frontier.cpu
#source /ccs/home/prlarsen/codes/HACC/env/bashrc.frontier.kokkos_hip
rm hydro_xray
rm *.o
make hydro_compile
make hydro_link



#!/bin/bash
source /ccs/home/prlarsen/codes/HACC/env/bashrc.frontier.cpu
#source /ccs/home/prlarsen/codes/HACC/env/bashrc.frontier.kokkos_hip
rm dens
rm *.o
make go_compile
make go_link



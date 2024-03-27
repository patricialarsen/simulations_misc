#!/bin/bash

#source /home/prlarsen/hacc_fresh/HACC/env/bashrc.polaris.kokkos_cuda
source /global/homes/p/plarsen/plarsen_git/HACC/env/bashrc.perlmutter.kokkos_cuda
#/home/prlarsen/hacc_fresh/HACC/env/bashrc.cooley.cpu
rm hydro_xray
rm *.o
make hydro_compile
make hydro_link



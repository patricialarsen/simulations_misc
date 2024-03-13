#!/bin/bash

source /home/prlarsen/hacc_fresh/HACC/env/bashrc.polaris.kokkos_cuda
#/home/prlarsen/hacc_fresh/HACC/env/bashrc.cooley.cpu
rm adiabatic_bhp
rm *.o
make bhp_compile
make bhp_link



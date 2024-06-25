#!/bin/bash

source /home/prlarsen/codes/HACC/env/bashrc.aurora.kokkos_sycl_thirdparty
#source /home/nfrontie/trunk/HACC/env/bashrc.aurora.kokkos_sycl

rm dens
rm *.o
make go_compile
make go_link



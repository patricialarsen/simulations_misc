#!/bin/bash

source /home/prlarsen/hacc_fresh/HACC/env/bashrc.cooley.cpu
rm hydro
rm *.o
make hydro_compile
make hydro_link



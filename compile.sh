#!/bin/bash

gfortran -g -ffpe-trap=zero,invalid,overflow,underflow differential_evolution.f90 graphf.f90 dijkstra.f90 utils_lj_fit.f90 lj_fit.f90 -o lj_fit.x
./lj_fit.x

#!/bin/sh

gfortran -c verlet_list.f90
gfortran -c forces_vlist_RaquelNEW.f90
gfortran -c program.f90

gfortran -o program verlet_list.o forces_vlist_RaquelNEW.o program.o

./program

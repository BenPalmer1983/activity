#!/bin/bash
gfortran -o activity.x src/kinds.f90 src/constants.f90 src/stringfunctions.f90 src/maths.f90 src/initialise.f90 src/input.f90 src/prep.f90 src/productionLoss.f90 src/output.f90 src/activity.f90 
sleep 1 
rm *.mod 
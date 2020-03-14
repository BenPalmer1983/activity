#!/bin/bash
#MOAB -l "nodes=1:ppn=2,walltime=0:05:00"
#MOAB -j oe
#MOAB -q bbtest
module load apps/openmpi/v1.6.3/gcc-tm-ib/v4.7.2
cd "/gpfs/bb/bxp912/activity/v1.0/examples"
mpirun -n 2 /gpfs/bb/bxp912/activity/v1.0/bin/activity.x -in Fe36MeVBB.in -out output.out

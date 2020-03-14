#!/bin/bash
#MOAB -l "nodes=4:ppn=8,walltime=2:00:00"
#MOAB -j oe
#MOAB -N pwscf_calc
#MOAB -A readmsd02
module load apps/openmpi/v1.6.3/gcc-tm-ib/v4.7.2
cd "$HOME/activity/v1.0/examples"
mpirun $HOME/activity/v1.0/bin/activity.x < Fe36MeVBB.in > output.out
#mpiexec -bynode -bysocket -bind-to-socket -display-map -report-bindings pw.x 
#!/bin/bash
# Load arguments to array
args=("$@")
# Change to current working directory (if not already done so)
cd ${0%/*}
# Store processor count
procCount="$(nproc 2>&1)"
#----------------------------------------------------------------------------------
# Prepare directory/file vars
workingDir=$(pwd)
defaultDir=$HOME"/activity"
dataFile=$workingDir"/data/activityData.tar.gz"
#----------------------------------------------------------------------------------
echo "=============================================================="
echo "Installing "$projectTitle
echo "University of Birmingham"
echo "=============================================================="
echo "Requirements: 500MB hard drive space"
echo "Working directory: "$workingDir
echo " "
echo "User Options:"
echo "------------------------------"
echo "Installation directory: (Default: $defaultDir)"
read dirIn
#----------------------------------------------------------------------------------
install=0
while [ $install -eq 0 ]; do
# Set default data directory
  if [[ $dirIn == "" || $dirIn == " " ]]
  then
    installationDir=$defaultDir
  else
    installationDir=$dirIn
  fi
# mk install dir
  mkdir -p $installationDir
# Check free space
  install=0
  result=$(df $installationDir)
# Read through result
  flag=0
  freeSpace=""
  for (( i=${#result}-1; i>0; i-- )); do
    if [ $flag == 0 ]
    then
      extract="${result:$i:2}"
      if [ "$extract" == "% " ]
      then
        flag=1
      fi
    elif [ $flag == 1 ]
    then
      extract="${result:$i:1}"
      if [ "$extract" == " " ]
      then
        flag=2
      fi
    elif [ $flag == 2 ]
    then
      extract="${result:$i:1}"
      if [ "$extract" != " " ]
      then
        flag=3
        freeSpace=$freeSpace$extract
      fi
    elif [ $flag == 3 ]
    then
      extract="${result:$i:1}"
      if [ "$extract" == " " ]
      then
        flag=4
      else
        freeSpace=$extract$freeSpace
      fi
    fi
  done
  if [ $freeSpace -gt 500000 ]
  then
    install=1
    freeSpaceM=$((freeSpace / 1000))
    echo "Free space: "$freeSpaceM"M.  Installing."
  else
    echo "Not enough space, please free space or choose another directory."
  fi
done
#----------------------------------------------------------------------------------
# set vars
binFile="activity.x"
binDir=$installationDir"/bin"
dataDir=$installationDir"/data"
examplesDir=$installationDir"/examples"
installLogDir=$installationDir"/installLog"
dataFileExtract=$dataDir"/activityData.tar.gz"
# Make directories
if [ -d "$installLogDir" ]; then
  rm -r $installLogDir
fi
mkdir -p $binDir
mkdir -p $dataDir
mkdir -p $examplesDir
mkdir -p $installLogDir
echo "Created installation directory "$installationDir
# Clear out install directory (if already existed)
#----------------------------------------------------------------------------------
# extract data
cp $dataFile $dataFileExtract
cd $dataDir
tar xzf activityData.tar.gz
#----------------------------------------------------------------------------------
# compile
cd $workingDir
modDir=$workingDir"/mod"

# make activity
echo $workingDir"/src"
cd "$workingDir"/src
mpif90 -O3 -g -Wall -mtune=native \
kinds.f90 \
logicalMod.f90 \
strings.f90 \
units.f90 \
constants.f90 \
general.f90 \
mpiSubs.f90 \
printMod.f90 \
arrayFunctions.f90 \
matrix.f90 \
linearAlgebra.f90 \
regression.f90 \
rng.f90 \
basicMaths.f90 \
solveFunctions.f90 \
laplaceTransforms.f90 \
calcFunctions.f90 \
activityFunctions.f90 \
msubs.f90 \
globals.f90 \
initialise.f90 \
input.f90 \
prep.f90 \
productionLoss.f90 \
output.f90 \
activity.f90 \
-J $modDir \
-o $binDir"/"$binFile

echo "Compiling activity.x:"
echo $commandLine
eval $commandLine
# add export to profile so user can run activity.x
exportLine="export PATH=\"\$PATH:"$binDir"\""
profileFile=$HOME"/.bash_profile"
touch $profileFile
if grep -q "$binDir" "$profileFile";
then
  echo $profileFile" already has the bin directory path"
else
  echo "Adding "$binDir" to "$profileFile
  echo $exportLine >> $profileFile
fi
source $profileFile
profileFile=$HOME"/.bashrc"
touch $profileFile
if grep -q "$binDir" "$profileFile";
then
  echo $profileFile" already has the bin directory path"
else
  echo "Adding "$binDir" to "$profileFile
  echo $exportLine >> $profileFile
fi
source $profileFile
#Make example file
cp $workingDir"/examples/Fe36MeV.exyz" $examplesDir"/Fe36MeV.exyz"
rm -f $examplesDir"/Fe36MeV.in"
touch $examplesDir"/Fe36MeV.in"
echo "#elements" >> $examplesDir"/Fe36MeV.in"
echo "Fe 100" >> $examplesDir"/Fe36MeV.in"
echo "#isotopes" >> $examplesDir"/Fe36MeV.in"
echo "\""$dataDir"/isotopes.txt\"" >> $examplesDir"/Fe36MeV.in"
echo "#decaymodes" >> $examplesDir"/Fe36MeV.in"
echo "\""$dataDir"/decaymodes.txt\"" >> $examplesDir"/Fe36MeV.in"
echo "#gammaenergies" >> $examplesDir"/Fe36MeV.in"
echo "\""$dataDir"/gammaenergies.txt\"" >> $examplesDir"/Fe36MeV.in"
echo "#xsfiles" >> $examplesDir"/Fe36MeV.in"
echo "\""$dataDir"/xs\"" >> $examplesDir"/Fe36MeV.in"
echo "#trajfile" >> $examplesDir"/Fe36MeV.in"
echo "\""$examplesDir"/Fe36MeV.exyz\"" >> $examplesDir"/Fe36MeV.in"
echo "#polyfitorder" >> $examplesDir"/Fe36MeV.in"
echo "5" >> $examplesDir"/Fe36MeV.in"
echo "#integrationgranularity" >> $examplesDir"/Fe36MeV.in"
echo "10" >> $examplesDir"/Fe36MeV.in"
echo "#beamflux" >> $examplesDir"/Fe36MeV.in"
echo "0.5 uA" >> $examplesDir"/Fe36MeV.in"
echo "#beamenergy" >> $examplesDir"/Fe36MeV.in"
echo "36 MeV" >> $examplesDir"/Fe36MeV.in"
echo "#beamduration" >> $examplesDir"/Fe36MeV.in"
echo "300 s" >> $examplesDir"/Fe36MeV.in"
echo "#beamarea" >> $examplesDir"/Fe36MeV.in"
echo "100 mm2" >> $examplesDir"/Fe36MeV.in"
echo "#amtime" >> $examplesDir"/Fe36MeV.in"
echo "260000 s" >> $examplesDir"/Fe36MeV.in"
echo "#timestep" >> $examplesDir"/Fe36MeV.in"
echo "1000 s" >> $examplesDir"/Fe36MeV.in"
echo "#projectile" >> $examplesDir"/Fe36MeV.in"
echo "1 1" >> $examplesDir"/Fe36MeV.in"
echo "#targetthickness" >> $examplesDir"/Fe36MeV.in"
echo "0.5 mm" >> $examplesDir"/Fe36MeV.in"
echo "#materialdensity" >> $examplesDir"/Fe36MeV.in"
echo "8000 kgm3" >> $examplesDir"/Fe36MeV.in"
echo "#vpi" >> $examplesDir"/Fe36MeV.in"
echo "60.2" >> $examplesDir"/Fe36MeV.in"
echo "#individualisotopeactivity" >> $examplesDir"/Fe36MeV.in"
echo "yes" >> $examplesDir"/Fe36MeV.in"
echo "#verboseterminal" >> $examplesDir"/Fe36MeV.in"
echo "yes" >> $examplesDir"/Fe36MeV.in"
echo "#targetdpa" >> $examplesDir"/Fe36MeV.in"
echo "0.0" >> $examplesDir"/Fe36MeV.in"
echo "#gammachartresolution" >> $examplesDir"/Fe36MeV.in"
echo "200" >> $examplesDir"/Fe36MeV.in"
echo "Installation complete"
echo "=============================================================="
echo "Install Complete"
echo "Running Example"
echo "=============================================================="
# Try running
exCommand="cd "$examplesDir"; mpirun -n 1 activity.x "$examplesDir"/Fe36MeV.in"
echo $exCommand
eval $exCommand
# Try running and save output
exCommand="cd "$examplesDir"; mpirun -n 1 activity.x "$examplesDir"/Fe36MeV.in > "$examplesDir"/Fe36MeV.out"
eval $exCommand





#

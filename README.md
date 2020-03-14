============================================================================
ACTIVITY code University of Birmingham
============================================================================

Requires:

Preprepared data files
Openmpi
Gfortran



Calculates the activity of a target irradiated by a proton beam.


Units
============================================================================
#beamflux                uA
#beamenergy              MeV
#beamduration            ms/s/m/hr/d
#beamarea                mm
#amtime                  ms/s/m/hr/d
#timestep                ms/s/m/hr/d          
#targetthickness         a/mm/cm/m
#materialdensity         kgm3/gcm3

Internally, the code uses uA for flux, barns for cross section, angstrom for depth, 
keV for energy, seconds for time, kgm-3 for density, number density in atoms per m3,
beam area in mm2.



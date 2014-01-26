Program activity

! Activity calculation code v0.02
! University of Birmingham
! Ben Palmer
!
! Internal units:
!   beam flux uA
!   cross sections barns
!   depth angstrom
!   energy keV
!   time s
!   number density atoms per m3
!   beam area mm2
!
! Requires exyz data file from SRIM
! Includes TALYS cross section data from TALYS endf data files
! Includes decay from Tuli + Sonzogni
! Isotope data from http://www.nist.gov/pml/data/comp.cfm


! Setup Modules
Use kinds				!data kinds
Use constants
Use stringfunctions		!string functions
Use maths
Use initialise			! input
Use input				! input
Use prep
Use productionLoss 
Use output				! input

Call runInitialise()
Call runInput()
Call runPrep()
Call runProductionLoss()
Call runOutput()




End
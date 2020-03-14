PROGRAM activity

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
  Use msubs
!  Use constants
!  Use stringfunctions    !string functions
!  Use maths
  Use globals
  Use initialise      ! input
  Use input        ! input
  Use prep
  Use productionLoss
  Use output        ! input

! Include MPI header
  Include 'mpif.h'

! Variables
  Integer(kind=StandardInteger) :: error
! Init MPI
  Call MPI_Init(error)
! Call subroutines
  Call initGlobals()
  Call runInitialise()
  Call runInput()
  Call runPrep()
  Call runProductionLoss()
  Call runOutput()

! Finalise MPI
  Call MPI_Finalize(error)

End Program activity

! --------------------------------------------------------------!
! MPI Subroutines module
! mpiSubsTypes, mpiSubs
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
! Makes configuration of atoms
! Calculates neighbour list for a collection of atoms
!
! ----------------------------------------
! Updated: 4th November 2015
! ----------------------------------------

Module mpiSubsTypes
! Setup Modules
  Use kinds
! Force declaration of all variables
  Implicit None
! Vars:  Module Parameters
!  Integer(kind=StandardInteger), Parameter ::
! Make private
  Private
! Public Variables and Parameters
!  Public ::
! Public derived types
  Public :: mpiObj

! Types
  Type :: mpiObj
    Integer(kind=StandardInteger) :: ID
    Integer(kind=StandardInteger) :: Count
  End Type mpiObj

End Module mpiSubsTypes


Module mpiSubs
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use mpi
  Use kinds
  Use strings
  Use general
  Use mpiSubsTypes
! Force declaration of all variables
  Implicit None
!---------------------------
! Variables
! Integer(kind=StandardInteger) :: i
! Make private
  Private
! Make Variables and Parameters Public
!  Public :: i
! Make Subroutines and Functions Public
  Public :: m_initMpi
  Public :: M_Loop

!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
!        Module Subroutines
!
! -----------------------------------------------

  Subroutine m_initMpi(mpiIn)
! Init the unit coords data type
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(mpiObj) :: mpiIn
! Vars:  Private
    Integer(kind=StandardInteger) :: error
! Init mpi obj
    Call MPI_Comm_size( MPI_COMM_WORLD,mpiIn%Count,error)
    Call MPI_Comm_rank(MPI_COMM_WORLD,mpiIn%ID,error)
  End Subroutine m_initMpi




! -----------------------------------------------
!        Module Functions
!
! -----------------------------------------------

  Function M_Loop(mpiIn, loop) Result (result)
! Checks if process should run on this loop
    Implicit None   ! Force declaration of all variables
! Vars:  In
    Type(mpiObj) :: mpiIn
    Integer(kind=StandardInteger) :: loop
! Vars:  Out
    Integer(kind=StandardInteger) :: loopID
    Logical :: result
! Assume false
    result = .false.
    loopID = Mod((loop-1),mpiIn%Count)
    If(mpiIn%ID.eq.loopID)Then
      result = .true.
    End If
  End Function M_Loop


End Module mpiSubs









!

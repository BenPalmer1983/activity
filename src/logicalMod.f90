Module logicalMod
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
! Force declaration of all variables
  Implicit None
! Make private
  Private
! Public
  Public :: FlipLogical


! Interfaces
!  Interface interfaceName
!    Module Procedure modNameA, modNameB, ...
!  End Interface interfaceName
!
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------

  Function FlipLogical(logicalIn) RESULT (logicalOut)
! If true return false, if false return true.
    Implicit None ! Force declaration of all variables
! In:      Declare variables
    Logical :: logicalIn
! Out:     Declare variables
    Logical :: logicalOut
! Private: Declare variables
    If(logicalIn)Then
      logicalOut = .false.
    Else
      logicalOut = .true.
    End If
  End Function FlipLogical


End Module logicalMod
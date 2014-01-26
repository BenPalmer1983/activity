Module kinds
  Implicit None
  Integer, Parameter :: sp = Selected_Real_Kind(6,37)    ! single real
  Integer, Parameter :: dp = Selected_Real_Kind(15,307)  ! double real
  Integer, Parameter :: qp = Selected_Real_Kind(33,4931) ! quadrupole real
  Integer, Parameter :: wp = dp                          ! working real
  Integer, Parameter :: ip = Selected_Int_Kind(12)       ! long integer
  Integer, Parameter :: ib = Selected_Int_Kind(32)       ! very long integer
End Module kinds
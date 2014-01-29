Module kinds
  Implicit None
  Integer, Parameter :: sp = Selected_Real_Kind(6,37)    ! single real
  Integer, Parameter :: dp = Selected_Real_Kind(15,307)  ! double real
  Integer, Parameter :: qp = Selected_Real_Kind(33,4931) ! quadrupole real
  Integer, Parameter :: wp = dp                          ! working real
  Integer, Parameter :: ip = Selected_Int_Kind(12)       ! long integer
  Integer, Parameter :: ib = Selected_Int_Kind(32)       ! very long integer
  
  Integer, Parameter :: SingleReal = Selected_Real_Kind(6,37)    	  ! single real, 6 decimal precision, exponent range 37		
  Integer, Parameter :: DoubleReal = Selected_Real_Kind(15,307)  	  ! double real, 15 decimal precision, exponent range 307		
  Integer, Parameter :: QuadrupoleReal = Selected_Real_Kind(33,4931)  ! quadrupole real
  Integer, Parameter :: SmallInteger = Selected_Int_Kind(4)           ! small integer
  Integer, Parameter :: StandardInteger = Selected_Int_Kind(8)        ! standard integer
  Integer, Parameter :: LongInteger = Selected_Int_Kind(12)           ! long integer
  Integer, Parameter :: VeryLongInteger = Selected_Int_Kind(32)       ! very long integer
End Module kinds
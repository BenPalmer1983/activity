Module kinds
  Implicit None
  Integer, Parameter :: SingleReal = Selected_Real_Kind(6,37)         ! single real, 6 decimal precision, exponent range 37
  Integer, Parameter :: DoubleReal = Selected_Real_Kind(15,307)       ! double real, 15 decimal precision, exponent range 307
  Integer, Parameter :: QuadrupoleReal = Selected_Real_Kind(33,4931)  ! quadrupole real
  Integer, Parameter :: TinyInteger = Selected_Int_Kind(1)            ! tiny integer         1 byte    -2^4 to 2^4-1                                           -16 to 16
  Integer, Parameter :: SmallInteger = Selected_Int_Kind(4)           ! small integer        4 bytes  -2^16 to 2^16-1                                       -65536 to 65535
  Integer, Parameter :: StandardInteger = Selected_Int_Kind(8)        ! standard integer     8 bytes  -2^31 to 2^31-1                                  -2147483648 to 2147483648
  Integer, Parameter :: LongInteger = Selected_Int_Kind(12)           ! long integer         16 bytes -2^63 to 2^63-1                         -9223372036854775808 to 9223372036854775807
  Integer, Parameter :: VeryLongInteger = Selected_Int_Kind(32)       ! very long integer    32 bytes -2^63 to 2^63-1     -170141183460469231731687303715884105728 to 170141183460469231731687303715884105727
End Module kinds

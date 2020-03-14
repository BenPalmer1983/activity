Module constants

! Setup Modules
  Use kinds

! force declaration of all variables
  Implicit None
! Useful parameters for source code
  Integer(kind=StandardInteger), Parameter :: maxFileRows = 10000000
! Physical and Mathematical constants
  Real(kind=DoubleReal), Parameter :: avogadrosConstant = 6.02214129D23
  Real(kind=DoubleReal), Parameter :: elementaryCharge = 1.60217656D-19
  Real(kind=DoubleReal), Parameter :: pi = 3.1415926535898D0
  Real(kind=DoubleReal), Parameter :: lnTwo = 0.6931471805599453094172321214581765680755001344D0
  Real(kind=DoubleReal), Parameter :: eulersNumber = 2.7182818284590452353602874713527D0
! Quadruple
  Real(kind=QuadrupoleReal), Parameter :: lnTwoQ = 0.6931471805599453094172321214581765680755001344D0

! Privacy of functions/subroutines/variables
  Private
! Useful parameters for source code
  Public :: maxFileRows
! Physical and Mathematical constants
  Public :: avogadrosConstant
  Public :: elementaryCharge
  Public :: pi
  Public :: lnTwo
  Public :: eulersNumber
  Public :: lnTwoQ

End Module constants
